
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





#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <gdal_utils.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <thread>
#include "timer.h"



// if true, input image will be cropped and subsampled in the first step at the same time
// as band selection and band rescaling, this adds an additional resammpling and hence might
// leas to slightly smoother results although this should not matter much in practice.
const bool ENABLE_EARLY_SUBSAMPLE = false;



//
struct arg_input {
    std::string output_path;
    std::vector<uint8_t> z;
    std::vector<uint32_t> x;
    std::vector<uint32_t> y;
    std::vector<uint16_t> b;
    std::vector<double> min;
    std::vector<double> max;
    std::string input_file;
    std::uint32_t m;
    bool enable_early_subsample;
    bool verbose;
} args;


typedef struct {
    double left;
    double right;
    double top;
    double bottom;
    std::string srs;
    uint32_t width;
    uint32_t height;
} space_view;


template<typename T>
struct bounds_2d {
    T left, bottom, top, right;
};

template<typename T>
struct coords_2d {
    T x, y;
};


const double EARTH_RADIUS_METERS = 6378137;
const uint16_t TILE_SIZE_PX = 256;

const double EARTH_CIRCUMFERENCE_METERS = 2 * M_PI * EARTH_RADIUS_METERS;
const double WEBMERCATOR_BOUNDS_LEFT = -EARTH_CIRCUMFERENCE_METERS / 2.0; // -20037508.342789244
const double WEBMERCATOR_BOUNDS_LOWER = -EARTH_CIRCUMFERENCE_METERS / 2.0; // -20037508.342789244
const double WEBMERCATOR_BOUNDS_RIGHT = EARTH_CIRCUMFERENCE_METERS / 2.0;  // 20037508.342789244
const double WEBMERCATOR_BOUNDS_UPPER = EARTH_CIRCUMFERENCE_METERS / 2.0;  // 20037508.342789244





/**
 * Returns tile boundaries in web mercator projection
 * @param tile_coords
 * @param zoom_level
 * @return
 */
bounds_2d<double> tile_bounds(coords_2d<uint32_t> tile_coords, uint8_t zoom_level) {
    bounds_2d<double> out;
    out.left = WEBMERCATOR_BOUNDS_LEFT + tile_coords.x * EARTH_CIRCUMFERENCE_METERS / std::pow(2, zoom_level);
    out.right = WEBMERCATOR_BOUNDS_LEFT + (tile_coords.x + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, zoom_level);
    out.top = WEBMERCATOR_BOUNDS_UPPER - tile_coords.y * EARTH_CIRCUMFERENCE_METERS / std::pow(2, zoom_level);
    out.bottom = WEBMERCATOR_BOUNDS_UPPER - (tile_coords.y + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, zoom_level);
    return out;
}


void affine(double a[6], double* iox, double* ioy, uint32_t n) {
    for (uint32_t i=0; i<n; ++i) {
        double x = a[0] + a[1] * iox[i] + a[2] * ioy[i];
        ioy[i] = a[3] + a[4] * iox[i] + a[5] * ioy[i];
        iox[i] = x;
    }
}


void affine_inv(double a[6], double* iox, double* ioy, uint32_t n) {
    double idet=1/(a[1]*a[5] - a[2]*a[4]);
    for (uint32_t i=0; i<n; ++i) {
        double x = idet *  ( a[5] * (iox[i] - a[0])   - a[2] * (ioy[i] - a[3])  );
        ioy[i] = idet *  ( - a[4] * (iox[i] - a[0])   + a[1] * (ioy[i] - a[3])  );
        iox[i] = x;
    }
}

int generate_tile(space_view sview, std::string filename) {


    sview.width = 256;
    sview.height = 256;



    // Step 1: Build in-memory VRT that has
    // - only selected input bands
    // - Byte data type and color scaling information

    std::string orig_vrt_path = args.input_file;

    GDALDataset *orig_vrt = (GDALDataset *) GDALOpenShared(orig_vrt_path.c_str(), GA_ReadOnly);
    if (orig_vrt == NULL) {
        std::cout << "ERROR: GDALOpen failed for original VRT." << std::endl;
        return 1;
    }


    std::string bandsel_vrt_path = "";
    CPLStringList translate_args(NULL);
    translate_args.AddString("-of");
    translate_args.AddString("VRT");

    for (uint8_t i = 0; i < args.b.size(); ++i) {
        translate_args.AddString("-b");
        translate_args.AddString(std::to_string(args.b[i]).c_str());
    }

    translate_args.AddString("-ot");
    translate_args.AddString("Byte");

    // TODO: try out gdal_translate with reduced resolution / cropped window (see below, not sure about performance benefit)
    if ( args.enable_early_subsample) {
        bounds_2d<double> dblcrop;
        dblcrop.left = sview.left - 0.5 * (sview.right - sview.left);
        dblcrop.right = sview.right + 0.5 * (sview.right - sview.left);
        dblcrop.top = sview.top + 0.5 * (sview.top - sview.bottom);
        dblcrop.bottom = sview.bottom - 0.5 * (sview.top - sview.bottom);

        translate_args.AddString("-projwin");
        translate_args.AddString(std::to_string(dblcrop.left).c_str());
        translate_args.AddString(std::to_string(dblcrop.top).c_str());
        translate_args.AddString(std::to_string(dblcrop.right).c_str());
        translate_args.AddString(std::to_string(dblcrop.bottom).c_str());
        translate_args.AddString("-projwin_srs");
        translate_args.AddString(sview.srs.c_str());

        translate_args.AddString("-outsize");
        translate_args.AddString(std::to_string(2*sview.width).c_str());
        translate_args.AddString(std::to_string(2*sview.height).c_str());
    }



    if (args.min.empty()) {
        translate_args.AddString("-scale");
    } else {
        for (uint16_t i = 0; i < args.b.size(); ++i) {
            translate_args.AddString((std::string("-scale_") + std::to_string(i + 1)).c_str());
            translate_args.AddString(std::to_string(args.min[i]).c_str());
            translate_args.AddString(std::to_string(args.max[i]).c_str());
        }
    }
    GDALTranslateOptions *trans_options = GDALTranslateOptionsNew(translate_args.List(), NULL);
    if (trans_options == NULL) {
        std::cout << "ERROR: Cannot create gdal_translate options." << std::endl;
        return 1;
    }

    GDALDatasetH bandsel_vrt = GDALTranslate(bandsel_vrt_path.c_str(), (GDALDatasetH) orig_vrt, trans_options, NULL);
    GDALTranslateOptionsFree(trans_options);








    // Step 2: Call gdalwarp to reproject and rescale the image based on tile parameters extent and zoom level
    // Result will be written to MEM dataset because PNG driver does not support Create (only CreateCopy).


    GDALDriverH mem_driver = GDALGetDriverByName("MEM"); // TODO: error handling
    GDALDatasetH warped_ds = GDALCreate(mem_driver, "", sview.width, sview.height, args.b.size() + 1, GDT_Byte,
                                        NULL); // additional alpha band
    GDALSetProjection(warped_ds, sview.srs.c_str());
    double affine[6];
    affine[0] = sview.left;
    affine[3] = sview.top;
    affine[1] = (sview.right - sview.left) / ((double) sview.width);
    affine[5] = (sview.bottom - sview.top) / ((double) sview.height);
    affine[2] = 0.0;
    affine[4] = 0.0;
    GDALSetGeoTransform(warped_ds, affine);

    CPLStringList warp_args(NULL);
    warp_args.AddString("-of");
    warp_args.AddString("MEM");
    warp_args.AddString("-t_srs");
    warp_args.AddString(sview.srs.c_str());

    warp_args.AddString("-te"); // xmin ymin xmax ymax
    warp_args.AddString(std::to_string(sview.left).c_str());
    warp_args.AddString(std::to_string(sview.bottom).c_str());
    warp_args.AddString(std::to_string(sview.right).c_str());
    warp_args.AddString(std::to_string(sview.top).c_str());

    warp_args.AddString("-te_srs");
    warp_args.AddString(sview.srs.c_str());

    warp_args.AddString("-ts");
    warp_args.AddString(std::to_string(sview.width).c_str());
    warp_args.AddString(std::to_string(sview.height).c_str());

    warp_args.AddString("-dstalpha");
    warp_args.AddString("-overwrite");
    warp_args.AddString("-multi");

    warp_args.AddString("-wm");
    warp_args.AddString(std::to_string(args.m).c_str());


    GDALWarpAppOptions *warp_opts = GDALWarpAppOptionsNew(warp_args.List(), NULL);
    if (warp_opts == NULL) {
        std::cout << "ERROR: Cannot create gdalwarp options." << std::endl;
        return 1;
    }
    GDALDatasetH inds = (GDALDatasetH) bandsel_vrt;

    //GDALWarp(NULL, warped_ds, 1,&inds, warp_opts,NULL);
    GDALWarp(NULL, warped_ds, 1, &inds, warp_opts, NULL);

    GDALWarpAppOptionsFree(warp_opts);




    // 3. Convert in-memory tile raster to PNG

    std::string png_path = filename;
    CPLStringList translate_png_args(NULL);
    translate_png_args.AddString("-of");
    translate_png_args.AddString("PNG"); // change to png

    GDALTranslateOptions *trans_png_options = GDALTranslateOptionsNew(translate_png_args.List(), NULL);
    if (trans_png_options == NULL) {
        std::cout << "ERROR: Cannot create gdal_translate options." << std::endl;
        return 1;
    }

    GDALTranslate(png_path.c_str(), warped_ds, trans_png_options, NULL);
    GDALTranslateOptionsFree(trans_png_options);




    // Step 4. clean up
    GDALClose(warped_ds);
}


int main(int ac, char *av[]) {




    // input args: input file, output directory, output format,  zoom levels, resampling,
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;


    po::options_description desc("Input arguments");
    desc.add_options()
            ("help,h", "print usage")
            ("output-dir", po::value<std::string>()->default_value(""),
             "output directory")
            ("zoom,z", po::value<std::string>(),
             "Zoom level(s)")
            ("tx,x", po::value<std::string>(),
             "tile x cooordinate(s)")
            ("ty,y", po::value<std::string>(),
             "tile y cooordinate(s)")
            ("bands,b", po::value<std::string>(),
             "Bands of the input dataset, either one or three must be selected")
            ("min", po::value<std::string>(),
             "minimum values that will be scaled to 0 for each band")
            ("max", po::value<std::string>(),
             "maximum values that will be scaled to 255 for each band")
            ("memory,m", po::value<uint32_t >()->default_value(256),
             "maximum memory in MiB for gdalwarp and GDAL cache")
            ("verbose,v", "verbose output")
            //("xx-ees", "enable early subsampling")
            ("input-file", po::value<std::string>());


    po::positional_options_description p;
    p.add("input-file", 1);
    p.add("output-dir", 1);
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

    args.enable_early_subsample = false;
//    if (vm.count("xx-ees")) {
//        args.enable_early_subsample = true;
//    }


    args.m = vm["memory"].as<uint32_t>();

    args.output_path = vm["output-dir"].as<std::string>();
    if (args.output_path.compare("") == 0) {
        fs::path p = fs::current_path();
        args.output_path = fs::absolute(p).string();
        std::cout << "INFO: using '" << args.output_path << "' as output directory" << std::endl;
    }
    fs::path out_path{args.output_path};
    if (!fs::exists(out_path)) {
        if (args.verbose) {
            std::cout << "Creating output directory '" << args.output_path << "." << std::endl;
        }
        // try to create
        if (!fs::create_directory(out_path)) {
            std::cout << "ERROR: failed to create output directory '" << args.output_path << "'." << std::endl;
        }
    }


    {
        std::istringstream str(vm["zoom"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint8_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint8_t) std::stoi(tokens[i]));
        }
        args.z = vec;
    }

    {
        std::istringstream str(vm["bands"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint16_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint16_t) std::stoi(tokens[i]));
        }
        args.b = vec;
    }

    {
        std::istringstream str(vm["tx"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint32_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint32_t) std::stoi(tokens[i]));
        }
        args.x = vec;
    }

    {
        std::istringstream str(vm["ty"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint32_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint32_t) std::stoi(tokens[i]));
        }
        args.y = vec;
    }

    if (vm.count("min") > 0) {
        std::istringstream str(vm["min"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<double> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((double) std::stod(tokens[i]));
        }
        args.min = vec;
    }

    if (vm.count("max") > 0) {
        std::istringstream str(vm["max"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<double> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((double) std::stod(tokens[i]));
        }
        args.max = vec;
    }


    if (!(args.z.size() == args.x.size() && args.x.size() == args.y.size())) {
        std::cout << "ERROR: input must have the same number of x,y,z tile coordinates." << std::endl;
    }
    if (args.z.size() < 1) {
        std::cout << "ERROR: at least one tile coordinate must be provided." << std::endl;
    }

    if (args.b.size() != 1 && args.b.size() != 3) {
        std::cout << "ERROR: either one or three bands must be selected." << std::endl;
    }
    if (args.min.size() != args.max.size()) {
        std::cout << "ERROR: expected exactly 0, 1, or 3 pairs of scaling parameters." << std::endl;
    }
    if (args.min.size() != args.b.size()) {
        if (args.min.size() == 1) { // use the same value for all bands
            args.min.push_back(args.min[0]);
            args.min.push_back(args.min[0]);
        } else if (args.min.size() != 0) {
            std::cout << "ERROR: expected 0, 1, or 3 scale (min and max) parameters." << std::endl;
        }
    }
    if (args.max.size() != args.b.size()) {
        if (args.max.size() == 1) { // use the same value for all bands
            args.max.push_back(args.max[0]);
            args.max.push_back(args.max[0]);
        } else if (args.max.size() != 0) {
            std::cout << "ERROR: expected 0, 1, or 3 scale parameters." << std::endl;
        }
    }


    if (!vm.count("input-file")) {
        std::cout << "ERROR: input file missing" << std::endl;
        return 1;
    }
    args.input_file = vm["input-file"].as<std::string>();




    GDALAllRegister();
    CPLSetConfigOption("GDAL_MAX_DATASET_POOL_SIZE", "400");
    GDALSetCacheMax(1024 * 1024 * args.m);
    CPLSetConfigOption("GDAL_PAM_ENABLED", "NO"); // avoid aux files for PNG tiles








    timer t, t1;
    t.start();
    for (uint16_t i = 0; i < args.z.size(); ++i) {
        t1.start();

        space_view sview;
        sview.srs = "EPSG:3857";

        coords_2d<uint32_t> tc;
        tc.x = args.x[i];
        tc.y = args.y[i];
        bounds_2d<double> ext = tile_bounds(tc, args.z[i]);


        sview.bottom = ext.bottom;
        sview.top = ext.top;
        sview.left = ext.left;
        sview.right = ext.right;
        sview.height = 256;
        sview.width = 256;

        std::string filename =
                std::to_string(args.z[i]) + "_" + std::to_string(args.x[i]) + "_" + std::to_string(args.y[i]) + ".png";
        generate_tile(sview, (fs::path{args.output_path} / fs::path{filename}).string());
        std::cout << "Tile (z/x/y) = (" << std::to_string(args.z[i]) << "," << args.x[i] << "," << args.y[i] << ") finished in " << t.time() << "s. " << std::endl;

    }
    std::cout << "DONE (" << t.time() << "s) ." << std::endl;
    return 0;
}