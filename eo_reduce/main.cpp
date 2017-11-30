
/*
 Copyright 2017 Marius Appel <marius.appel@uni-muenster.de>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/





/*
 *
 * eo_reduce -o /home/marius/Desktop/test_reduce.tif /home/marius/github/eo_webtools/sentinel_webview/data/S2A_MSIL1C_20161213T093402_N0204_R136_T34UFD_20161213T093403.vrt /home/marius/github/eo_webtools/sentinel_webview/data/S2A_MSIL1C_20170411T081601_N0204_R121_T34KGD_20170411T083418.vrt /home/marius/github/eo_webtools/sentinel_webview/data/S2A_MSIL1C_20171002T094031_N0205_R036_T34UFD_20171002T094035.vrt
 * eo_reduce -f "mean" --t_win="799500 809500 789000 799000" --t_size="334 334" --t_srs="EPSG:32636" -o /home/marius/Desktop/test_reduce_landsat.tif /home/marius/github/eo_webtools/eo_reduce/data/Landsat/
 *
 *
 *
 *  v.left = 799500.000;
    v.right = 809500.000;
    v.top = 799000.000;
    v.bottom = 789000.000;
    v.width = 334;
    v.height = 334;
    v.srs = "EPSG:32636";

 *
 */



/*
 * TODO:
 * - add resampling argument
 * - implement band selection
 * - add view argument
 * - add memorylimit argument
 */




#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <gdal_utils.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <thread>
#include <queue>
#include "timer.h"



#define BLOCK_SIZE 256
#define MEM_HARD_LIMIT_MB 512

typedef struct {
    double left;
    double right;
    double top;
    double bottom;
    std::string srs;
    uint32_t width;
    uint32_t height;
} space_view;


struct arg_input {
    std::string output_path; // output file (TIF + VRT?)
    std::vector<std::string> input_image;
    bool verbose;
    std::vector<uint16_t> b; // bands
    uint16_t m;
    space_view sview;
    double (*func)(double*, uint32_t, uint32_t);
} args;




template<typename T>
struct bounds_2d {
    T left, bottom, top, right;
};

template<typename T>
struct coords_2d {
    T x, y;
};







double median(double *buf, uint32_t delta, uint32_t n) {

    if (n==0) return nan("");

    std::vector<double> v;
    for (uint32_t i=1; i<n; ++i)
        v.push_back(buf[i * delta]);

    std::sort(v.begin(), v.end());
    if (v.size() % 2 == 0) {
        return (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2.0;
    }
    else return v[v.size() / 2];

//    std::priority_queue<double, std::priority_queue<double>::container_type, std::less<double> > left;
//    std::priority_queue<double, std::priority_queue<double>::container_type, std::greater<double> > right;
//
//    if (n==0) return nan("");
//    left.push(buf[0]);
//
//    for (uint32_t i=1; i<n; ++i) {
//        double el = buf[i * delta];
//
//        if (left.size() > right.size()) {
//            double l=left.top();
//            if (el < l) {
//                left.pop();
//                right.push(l);
//                left.push(el);
//            }
//            else right.push(el);
//        }
//        else {
//            double r = right.top();
//            if (el > r) {
//                right.pop();
//                left.push(r);
//                right.push(el);
//            }
//            else {
//                left.push(el);
//            }
//        }
//        assert(left.size() + right.size() == i+1);
//    }
//    if (right.size() == left.size()) return (right.top() + left.top()) / 2.0;
//    return left.top();
}





inline double mean(double *buf, uint32_t delta, uint32_t n) {
    if (n==0) return nan("");
    uint32_t count=0;
    double sum = 0.0;
    for (uint32_t i=0; i<n; ++i) {
        if (std::isfinite(buf[i * delta])) {
            ++count;
            sum += buf[i * delta];
        }
    }
    return sum / count;
}






double min(double *buf, uint32_t delta, uint32_t n) {
    if (n==0) return nan("");
    double min = DBL_MAX;
    for (uint32_t i=0; i<n; ++i) {
        if (std::isfinite(buf[i * delta])) {
            if (buf[i * delta] < min) min = buf[i * delta];
        }
    }
    return std::isfinite(min)? min : nan("");
}


double max(double *buf, uint32_t delta, uint32_t n) {
    if (n==0) return nan("");
    double max = -DBL_MAX;
    for (uint32_t i=0; i<n; ++i) {
        if (std::isfinite(buf[i * delta])) {
            if (buf[i * delta] > max) max = buf[i * delta];
        }
    }
    return std::isfinite(max)? max : nan("");
}





bool intersects(bounds_2d<double> a, bounds_2d<double> b) {
    return (
        a.right >= b.left &&
        a.left <= b.right &&
        a.top >= b.bottom &&
        a.bottom <= b.top);
}



inline uint32_t buf_index(uint32_t ii, uint16_t ib, uint32_t ix, uint32_t iy,
                          uint32_t ni, uint16_t nb, uint32_t nx, uint32_t ny) {
    //return ii * (nb * nx * ny) + ib * (nx * ny) + iy * (nx) + ix * 1; // assuming order x fastest, y 2nd, band 3rd, image 4th
    return ii * (nb * nx * ny) + iy * (nx * nb) + ix * (nb) + ib * 1; // assuming order band fastest, x 2nd, y 3rd, image 4th
}

// returns the gap in the buffer between pixel values from different images but equal x, y, and band
inline uint32_t buf_gap_image(uint32_t ni, uint16_t nb, uint32_t nx, uint32_t ny) {
    // this should return ni;
    return buf_index(1,0,0,0,ni,nb,nx,ny) - buf_index(0,0,0,0,ni,nb,nx,ny);
}

// returns the gap in the buffer between pixel values from different bands but equal x, y, and image
inline uint32_t buf_gap_band(uint32_t ni, uint16_t nb, uint32_t nx, uint32_t ny) {
    return buf_index(0,1,0,0,ni,nb,nx,ny) - buf_index(0,0,0,0,ni,nb,nx,ny);
}

inline uint32_t buf_gap_x(uint32_t ni, uint16_t nb, uint32_t nx, uint32_t ny) {
    return buf_index(0,0,1,0,ni,nb,nx,ny) - buf_index(0,0,0,0,ni,nb,nx,ny);
}

inline uint32_t buf_gap_y(uint32_t ni, uint16_t nb, uint32_t nx, uint32_t ny) {
    return buf_index(0,0,0,1,ni,nb,nx,ny) - buf_index(0,0,0,0,ni,nb,nx,ny);
}



int run() {



    // 1. Extract spatial extent of all input datasets

    // 2. Create empty output dataset (GTiff?)

    // 3. Over all blocks of the output dataset (e.g. 256x256),

    // 3.1 gather relevant scenes

    // 3.2 Call gdalwarp and read data to MEM datasets

    // 3.3 call RasterIO for all scenes

    // 3.4. Apply reduce function over returned arrays over all bands

    // 3.5. write output block of input image

    // 4. [Create overviews]

    // 5. [Create VRT]

    // further ideas:
    // - cache chunks of input data
    // - hold pool for open input datasets to prevent out of open file limits from OS


    // ASSUMPTIONS
    // a) All images have the same number of bands and bands from different datasets with the same band number have the same data type
    // b) Output datatype is double for all bands


    OGRSpatialReference srs_out;
    srs_out.SetFromUserInput(args.sview.srs.c_str());
    double affine_out[6];
    affine_out[0] = args.sview.left;
    affine_out[3] = args.sview.top;
    affine_out[1] = (args.sview.right - args.sview.left) / (double)args.sview.width;
    affine_out[5] = -(args.sview.top - args.sview.bottom) / (double)args.sview.height;
    affine_out[2] = 0;
    affine_out[4] = 0;

    GDALDriverH mem_driver = GDALGetDriverByName("MEM");
    GDALDriverH gtiff_driver = GDALGetDriverByName("GTiff");


    uint16_t n_bands = 0;
    int* band_list = NULL;

    std::vector<bounds_2d<double>> extent;
    for (uint32_t i=0; i<args.input_image.size(); ++i) {
        GDALDataset *gdal_in = (GDALDataset *) GDALOpenShared(args.input_image[i].c_str(), GA_ReadOnly);
        if (gdal_in == NULL) {
            std::cout << "ERROR: GDALOpen failed for '" << args.input_image[i] << "'" <<  std::endl;
            return 1;
        }

        if (i==0) {
            n_bands = gdal_in->GetRasterCount(); // ASSUMPTION a)
            band_list = (int*)calloc(n_bands, sizeof(int));
            for (uint16_t ib=0; ib<n_bands; ++ib) band_list[ib] = ib+1;
        }


        OGRSpatialReference srs_in(gdal_in->GetProjectionRef()); // Does this work?
        double affine_in[6];
        gdal_in->GetGeoTransform(affine_in);

        OGRCoordinateTransformation *coord_transform = OGRCreateCoordinateTransformation( &srs_in, &srs_out );

        double x[4] = {affine_in[0],
                          affine_in[0],
                          affine_in[0] + affine_in[1] * gdal_in->GetRasterXSize() + affine_in[2] * gdal_in->GetRasterYSize(),
                          affine_in[0] + affine_in[1] * gdal_in->GetRasterXSize() + affine_in[2] * gdal_in->GetRasterYSize()};
        double y[4] = {affine_in[3],
                          affine_in[3] + affine_in[4] * gdal_in->GetRasterXSize() + affine_in[5] * gdal_in->GetRasterYSize(),
                          affine_in[3],
                          affine_in[3] + affine_in[4] * gdal_in->GetRasterXSize() + affine_in[5] * gdal_in->GetRasterYSize()};



        if( coord_transform == NULL || !coord_transform->Transform( 4, x, y )) {
            std::cout << "ERROR: coordinate transformation failed" << std::endl;
            return 1;
        }

        double xmin = DBL_MAX;
        double ymin = DBL_MAX;
        double xmax = -DBL_MAX;
        double ymax = -DBL_MAX;
        for (uint8_t k=0; k<4; ++k) {
            if (x[k] < xmin ) xmin = x[k];
            if (y[k] < ymin ) ymin = y[k];
            if (x[k] > xmax ) xmax = x[k];
            if (y[k] > ymax ) ymax = y[k];
        }

        bounds_2d<double> in_extent;
        in_extent.left = xmin;
        in_extent.right = xmax;
        in_extent.top = ymax;
        in_extent.bottom = ymin;

        extent.push_back(in_extent);

        // TODO: derive spatial index?
        // STRTree envelope -> int (index in vector)

        GDALClose((GDALDatasetH)gdal_in);
    }




    CPLStringList out_co(NULL);
    out_co.AddNameValue("TILED","YES");
    out_co.AddNameValue("BLOCKXSIZE","256");
    out_co.AddNameValue("BLOCKYSIZE","256");
    out_co.AddNameValue("COMPRESS","DEFLATE");


    GDALDatasetH gdal_out = GDALCreate(gtiff_driver, args.output_path.c_str(), args.sview.width, args.sview.height, n_bands, GDT_Float64, out_co.List());
    char *wkt_out;
    srs_out.exportToWkt(&wkt_out);
    GDALSetProjection(gdal_out, wkt_out);
    GDALSetGeoTransform(gdal_out,affine_out);





    uint32_t n_blocks_x = (uint32_t)std::ceil((double)args.sview.width / (double)BLOCK_SIZE);
    uint32_t n_blocks_y = (uint32_t)std::ceil((double)args.sview.height / (double)BLOCK_SIZE);

    for (uint32_t ix=0; ix < n_blocks_x; ++ix) {
        for (uint32_t iy=0; iy < n_blocks_y; ++iy) {

            // Get extent of block
            bounds_2d<double> block_extent;
            block_extent.bottom = affine_out[3] + affine_out[5] * fmin((iy+1)*BLOCK_SIZE, args.sview.height);
            block_extent.top    = affine_out[3] + affine_out[5] * (iy)*BLOCK_SIZE;
            block_extent.right = affine_out[0] + affine_out[1] * fmin((ix+1)*BLOCK_SIZE, args.sview.width) ;
            block_extent.left = affine_out[0] + affine_out[1] * (ix)*BLOCK_SIZE;

            bounds_2d<uint32_t>block_extent_int; // direction of axes is not important here
            block_extent_int.left =  ix * BLOCK_SIZE;
            block_extent_int.top  =  iy * BLOCK_SIZE;
            block_extent_int.right =  fmin(block_extent_int.left + BLOCK_SIZE, args.sview.width);
            block_extent_int.bottom =  fmin(block_extent_int.top + BLOCK_SIZE, args.sview.height);

            // this is the actual width and height of the current block, which potentially has been cropped to
            // the extent of the output image
            uint32_t block_width =std::abs(block_extent_int.right - block_extent_int.left);
            uint32_t block_height = std::abs(block_extent_int.bottom - block_extent_int.top);


            // Gather all scenes which intersect with the current block
            // TODO: use spatial index to make filtering more efficient
            std::vector<uint32_t> image_ids;
            for (uint32_t ii=0; ii<args.input_image.size(); ++ii) {
                if (intersects(block_extent,extent[ii])) {
                    image_ids.push_back(ii);
                   // std::cout << "Block (" << ix << "," << iy << ") intersects with image " << ii << "."  << std::endl;
                }
               // else std::cout << "Block (" << ix << "," << iy  << ") DOES NOT intersect with image " << ii << "."  << std::endl;
                //print_extent(block_extent);
               // print_extent(extent[ii]);
            }



            // allocate memory to read data from all bands
            uint64_t buf_size_in = block_width*block_height*image_ids.size() * n_bands * sizeof(double);
            if (buf_size_in / (1024 * 1024) > MEM_HARD_LIMIT_MB) {
                std::cout << "ERROR: reached memory limit of " << MEM_HARD_LIMIT_MB << " MiB, please choose smaller block sizes." <<  std::endl;
                return 1;
            }
            double* buf_in = (double*)calloc(block_width*block_height*image_ids.size() * n_bands, sizeof(double));

            // TODO not needed?
            for (uint32_t j=0; j<block_width*block_height*image_ids.size() * n_bands; ++j) {
                buf_in[j] = nan("");
            }

            // Read data for block from all scenes
            for (uint32_t ii=0; ii<image_ids.size(); ++ii) {

                GDALDataset *gdal_in = (GDALDataset *) GDALOpenShared(args.input_image[image_ids[ii]].c_str(), GA_ReadOnly);
                if (gdal_in == NULL) {
                    std::cout << "ERROR: GDALOpen failed for '" << args.input_image[image_ids[ii]] << "'" <<  std::endl;
                    return 1;
                }


                // Call GDALWarp
                GDALDatasetH gdal_warped = GDALCreate(mem_driver, "", block_width, block_height, n_bands, GDT_Float64, NULL);
                //GDALDatasetH gdal_warped = GDALCreate(gtiff_driver, ("/home/marius/Desktop/" + std::to_string(ii) + ".tif" ).c_str(), block_width, block_height, n_bands, GDT_Float64, NULL);
                for (uint16_t ib=0; ib<n_bands; ++ib) {
                    GDALRasterBandH b = GDALGetRasterBand(gdal_warped, ib+1);
                    GDALSetRasterNoDataValue(b, nan(""));
                    //GDALFillRaster(b,nan(""), 0);
                }

                GDALSetProjection(gdal_warped, wkt_out);
                double affine[6];
                affine[0] = block_extent.left;
                affine[3] = block_extent.top;
                affine[1] = (block_extent.right - block_extent.left) / ((double) block_width);
                affine[5] = (block_extent.bottom - block_extent.top) / ((double) block_height);
                affine[2] = 0.0;
                affine[4] = 0.0;
                GDALSetGeoTransform(gdal_warped, affine);

                CPLStringList warp_args(NULL);
                warp_args.AddString("-of");
                warp_args.AddString("MEM");
                warp_args.AddString("-t_srs");
                warp_args.AddString(args.sview.srs.c_str());

                warp_args.AddString("-dstnodata");
                warp_args.AddString("nafn");

                warp_args.AddString("-te"); // xmin ymin xmax ymax
                warp_args.AddString(std::to_string(block_extent.left).c_str());
                warp_args.AddString(std::to_string(block_extent.bottom).c_str());
                warp_args.AddString(std::to_string(block_extent.right).c_str());
                warp_args.AddString(std::to_string(block_extent.top).c_str());

                warp_args.AddString("-te_srs");
                warp_args.AddString(args.sview.srs.c_str());

                warp_args.AddString("-ts");
                warp_args.AddString(std::to_string(block_width).c_str());
                warp_args.AddString(std::to_string(block_height).c_str());



                GDALWarpAppOptions *warp_opts = GDALWarpAppOptionsNew(warp_args.List(), NULL);
                if (warp_opts == NULL) {
                    std::cout << "ERROR: Cannot create gdalwarp options." << std::endl;
                    return 1;
                }
                GDALDatasetH gdal_warpin = (GDALDatasetH) gdal_in;
                GDALWarp(NULL, gdal_warped, 1, &gdal_warpin, warp_opts, NULL);
                GDALWarpAppOptionsFree(warp_opts);

                GDALDatasetRasterIOEx(gdal_warped,GF_Read, 0, 0, block_width, block_height, &buf_in[ buf_index(ii, 0,0,0,image_ids.size(), n_bands, block_width, block_height)],
                                      block_width, block_height, GDT_Float64, n_bands, NULL, 0, 0, 0, NULL);

                GDALClose(gdal_warped);
                GDALClose((GDALDatasetH)gdal_in);
            }





            uint32_t delta_i = buf_gap_image(image_ids.size(), n_bands,block_width, block_height);
            uint32_t delta_b = buf_gap_band(image_ids.size(), n_bands,block_width, block_height);
            uint32_t delta_x = buf_gap_x(image_ids.size(), n_bands,block_width, block_height);
            uint32_t delta_y = buf_gap_y(image_ids.size(), n_bands,block_width, block_height);



            // Do the analysis
            // uint64_t buf_size_out= BLOCK_SIZE*BLOCK_SIZE * n_bands * sizeof(double);
            double* buf_out = (double*)calloc(block_width*block_height * n_bands,sizeof(double));


            // test: fill buf out with nan
            for (uint32_t j=0; j<block_width*block_height * n_bands; ++j) {
                buf_out[j] = nan("");
            }

            if (!image_ids.empty()) {
                for (uint16_t in_ib = 0; in_ib < n_bands; ++in_ib) {
                    for (uint32_t in_ix = 0; in_ix < block_width; ++in_ix) {
                        for (uint32_t in_iy = 0; in_iy < block_height; ++in_iy) {
                            uint32_t i = buf_index(0, in_ib, in_ix, in_iy, image_ids.size(), n_bands, block_width,
                                                   block_height);

                            double val = args.func(&buf_in[i], delta_i, image_ids.size());


                           // std::cout << i << " in [" << 0 << "," << block_width * block_height * image_ids.size() * n_bands << "]" << std::endl;
                            if (!std::isfinite(val)) val = nan("");
                            buf_out[buf_index(0, in_ib, in_ix, in_iy, 1, n_bands, block_width, block_height)] = val;

                        }
                    }
                    // maybe write single band already here
                }
            }


            GDALDatasetRasterIOEx(gdal_out, GF_Write, block_extent_int.left, block_extent_int.top, block_width,block_height, buf_out, block_width, block_height, GDT_Float64, n_bands, NULL, 0, 0, 0, NULL);



            free(buf_in);
            free(buf_out);
           // if (band_list != NULL) free(band_list);
        }
    }



    for (uint16_t ib=0; ib<n_bands; ++ib) {
        GDALRasterBandH b = GDALGetRasterBand(gdal_out, ib+1);
        GDALSetRasterNoDataValue(b,nan(""));
    }

    GDALClose((GDALDatasetH)gdal_out);

}






int main(int ac, char *av[]) {


    // input args: input file, output directory, output format,  zoom levels, resampling,
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;


    po::options_description desc("Input arguments");
    desc.add_options()
            ("help,h", "print usage")
            ("output-file,o", po::value<std::string>(),
             "path of output file")
            ("bands,b", po::value<std::string>(),
             "Bands of the input dataset, defaults to all available bands")
            ("t_win", po::value<std::string>(), "coordinates of the target window expressed in the target reference system  with order (xmin, xmax, ymin, ymax).")
            ("t_size", po::value<std::string>(), "size of the target image in pixels with order (width, height).")
            ("t_srs", po::value<std::string>(), "target reference system as GDAL readable string.")
            ("func,f", po::value<std::string>(), "reducer function, currently supported are 'mean, 'min', 'max', 'median")
            ("memory,m", po::value<uint16_t>()->default_value(256), "memory in MB that is available for caching")
            ("verbose,v", "verbose output")
            ("input-file", po::value<std::vector<std::string>>());


    po::positional_options_description p;
    p.add("input-file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).
            options(desc).positional(p).run(), vm);
    po::notify(vm);


    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    args.verbose = false;
    if (vm.count("verbose")) {
        args.verbose = true;
    }


    args.m = vm["memory"].as<uint16_t>();



    if (!vm.count("output-file"))
    {
        std::cout << "ERROR: missing argument output file, specify with -o \"/PATH/TO/FILE\""  << std::endl;
        return 1;
    }
    args.output_path = vm["output-file"].as<std::string>();
    if (fs::exists( args.output_path)) {
        std::cout << "INFO: already existing output file will be overwritten." << std::endl;
    }



    if (vm.count("bands"))
    {
        std::istringstream str(vm["bands"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint16_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint16_t) std::stoi(tokens[i]));
        }
        args.b = vec;
    }
    else {
        // TODO fill bands
    }

    if (!vm.count("t_win"))
    {
        std::cout << "ERROR: no spatial target window given." << std::endl;
        return 1;
    }
    else {
        std::istringstream str(vm["t_win"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<double> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((double) std::stod(tokens[i]));
        }
        if (vec.size() != 4) {
            std::cout << "ERROR: invalid --t_win option provided." << std::endl;
            return 1;
        }
        args.sview.left = vec[0];
        args.sview.right = vec[1];
        args.sview.bottom = vec[2];
        args.sview.top = vec[3];
    }

    if (!vm.count("t_size"))
    {
        std::cout << "ERROR: no spatial target size given." << std::endl;
        return 1;
    }
    else {
        std::istringstream str(vm["t_size"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint32_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint32_t) std::stoi(tokens[i]));
        }
        if (vec.size() != 2) {
            std::cout << "ERROR: invalid --t_size option provided." << std::endl;
            return 1;
        }
        args.sview.width = vec[0];
        args.sview.height = vec[1];
    }

    if (!vm.count("t_srs"))
    {
        std::cout << "ERROR: no target spatial reference system given." << std::endl;
        return 1;
    }
    else {
        args.sview.srs = vm["t_srs"].as<std::string>();
    }


    if (!vm.count("func"))
    {
        std::cout << "ERROR: no reducer function provided." << std::endl;
        return 1;
    }
    else {
        std::string func = vm["func"].as<std::string>();
        std::transform(func.begin(),func.end(),func.begin(),::tolower);
        if(func.compare("mean") == 0 || func.compare("avg") == 0) args.func = &mean;
        else if(func.compare("median") == 0) args.func = &median;
        else if(func.compare("min") == 0) args.func = &min;
        else if(func.compare("max") == 0) args.func = &max;
        else {
            std::cout << "ERROR: invalid reducer function given." << std::endl;
            return 1;
        }
    }




    if (!vm.count("input-file")) {
        std::cout << "ERROR: missing input files or input directory" << std::endl;
        return 1;
    }

    args.input_image = vm["input-file"].as<std::vector<std::string>>();

    if (vm.count("input-file") == 1) {
        fs::path p(args.input_image[0]);

        if (fs::is_directory(p)) {
            args.input_image.clear();
            for(auto& entry : boost::make_iterator_range(fs::directory_iterator(p), {})) {
                args.input_image.push_back(entry.path().string());
            }
            std::cout << "INFO: Added " << args.input_image.size() << " images to input collection." << std::endl;
        }
    }



    GDALAllRegister();
    CPLSetConfigOption("GDAL_MAX_DATASET_POOL_SIZE", "400");
    //CPLSetConfigOption("GDAL_RASTERIO_RESAMPLING", "NEAR"); // does this work? if yes, even for RasterIO (non Ex) and gdalwarp?
    GDALSetCacheMax(1024 * 1024 * args.m);





    timer t;
    t.start();

    run();

    std::cout << "DONE (" << t.time() << "s) ." << std::endl;




    return 0;
}