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


struct arg_input {
    std::string output_file;
    std::string input_file;
    bool cog;
    bool verbose;
    std::vector<uint16_t> overview_levels;
} args;





int create_vrt() {



    GDALDriver *vrt_driver = (GDALDriver*) GDALGetDriverByName("VRT");
    if (vrt_driver == NULL) {
        std::cout << "ERROR: cannot find GDAL driver for VRT." << std::endl;
        return 1;
    }


    GDALDataset *in_dataset = (GDALDataset *)GDALOpenShared(args.input_file.c_str(), GA_ReadOnly);
    if (in_dataset == NULL) {
        std::cout << "ERROR: GDALOpen failed" << std::endl;
        return 1;
    }



    std::map<std::string, std::map<std::string, std::string> > scene_sds;



    char **sd_md = in_dataset->GetMetadata("SUBDATASETS");
    uint8_t n_sd = CSLCount(sd_md) / 2; // NAME + DESCR per SUBDATASET
    std::vector<std::string> sd;
    for (uint8_t i=0; i<n_sd; ++i) {
        std::string key = "SUBDATASET_" + std::to_string(i+1) + "_NAME";
        std::string val =  CSLFetchNameValue(sd_md, key.c_str());


        int32_t p=0;
        std::vector<std::string> v;
        while (p != std::string::npos) {
            int32_t p1 = val.find(":",p);
            if (p1 != std::string::npos) {
                v.push_back(val.substr(p,p1-p));
                p = p1 + 1;
            }
            else {
                v.push_back(val.substr(p,val.length()-p));
                break;
            }
        }

        assert(v.size() == 4);
        // v[0] is e.g. SENTINEL2_L1C
        // v[1] is path
        // v[2] is resolution or TCI
        // v[3] is EPSG_XXXXX

        if (scene_sds.find(v[3]) == scene_sds.end()) {
            std::map<std::string, std::string> dummy;
            scene_sds.insert(std::pair<std::string, std::map<std::string, std::string>>(v[3], dummy));
        }
        scene_sds[v[3]].insert(std::pair<std::string, std::string>(v[2],val));


        //sd.push_back();
        //std::cout << sd[i] << std::endl;
    }

    if (scene_sds.size() > 1) {
        std::cout << "WARNING: the scene contains images from different UTM zones, separate VRT datasets will be created.";
    }




    std::map<std::string, std::map<std::string, std::string> >::iterator it = scene_sds.begin();

    uint16_t is = 0;
    while (it != scene_sds.end()) {
        // TODO: assert that 10m, 20m, and 60m subdatasets are available

        GDALDataset *in_sd_10m = (GDALDataset *) GDALOpenShared(it->second["10m"].c_str(), GA_ReadOnly);
        GDALDataset *in_sd_20m = (GDALDataset *) GDALOpenShared(it->second["20m"].c_str(), GA_ReadOnly);
        GDALDataset *in_sd_60m = (GDALDataset *) GDALOpenShared(it->second["60m"].c_str(), GA_ReadOnly);
        if (in_sd_10m == NULL) {
            std::cout << "ERROR: GDALOpen failed for subdataset " << it->second["10m"].c_str() << std::endl;
            return 1;
        }
        if (in_sd_20m == NULL) {
            std::cout << "ERROR: GDALOpen failed for subdataset " << it->second["20m"].c_str() << std::endl;
            return 1;
        }
        if (in_sd_60m == NULL) {
            std::cout << "ERROR: GDALOpen failed for subdataset " << it->second["60m"].c_str() << std::endl;
            return 1;
        }




        if (args.cog) {
            if (args.verbose) {
                std::cout << "... converting scene to cloud-optimized GeoTIFF (" << is + 1 << "/" << scene_sds.size() << ")" << std::endl;
            }
            // -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE
            CPLStringList translate_args(NULL);
            translate_args.AddString("-of");
            translate_args.AddString("GTiff");

            translate_args.AddString("-co");
            translate_args.AddString("TILED=YES");

            translate_args.AddString("-co");
            translate_args.AddString("COPY_SRC_OVERVIEWS=YES");

            translate_args.AddString("-co");
            translate_args.AddString("COMPRESS=DEFLATE");

            translate_args.AddString("-co");
            translate_args.AddString("ZLEVEL=6");



            GDALTranslateOptions *trans_options = GDALTranslateOptionsNew(translate_args.List(), NULL);
            if (trans_options == NULL) {
                std::cout << "ERROR: Cannot create gdal_translate options." << std::endl;
                return 1;
            }

            // TODO: where to store COGs?

            GDALDatasetH tmp;
            boost::filesystem::path dir = boost::filesystem::path{args.input_file}.parent_path();
            tmp = GDALTranslate((dir / ( "10m" + it->first + ".tif")).string().c_str(), (GDALDatasetH) in_sd_10m, trans_options, NULL);
            GDALClose(in_sd_10m);
            in_sd_10m = (GDALDataset*)tmp;
            tmp = GDALTranslate((dir / ( "20m" + it->first + ".tif")).string().c_str(), (GDALDatasetH) in_sd_20m, trans_options, NULL);
            GDALClose(in_sd_20m);
            in_sd_20m = (GDALDataset*)tmp;
            tmp = GDALTranslate((dir / ( "60m" + it->first + ".tif")).string().c_str(), (GDALDatasetH) in_sd_60m, trans_options, NULL);
            GDALClose(in_sd_60m);
            in_sd_60m = (GDALDataset*)tmp;


            GDALTranslateOptionsFree(trans_options);

        }






        // TODO: assert that subdatasets have same extent
        std::string wkt_10m = in_sd_10m->GetProjectionRef();
        std::string wkt_20m = in_sd_20m->GetProjectionRef();
        std::string wkt_60m = in_sd_60m->GetProjectionRef();

        assert(wkt_10m.compare(wkt_20m) == 0 && wkt_20m.compare(wkt_60m) == 0);

        // Create target dataset
        std::string outfile = (scene_sds.size() == 1)? args.output_file : (args.output_file + "_" + it->first);
        GDALDataset *out_vrt = vrt_driver->Create(outfile.c_str(), in_sd_10m->GetRasterXSize(),in_sd_10m->GetRasterYSize(), 0, GDT_Unknown,NULL );


        GDALRasterBand* in_bands[13];
        in_bands[0] = in_sd_60m->GetRasterBand(1);
        in_bands[1] = in_sd_10m->GetRasterBand(1);
        in_bands[2] = in_sd_10m->GetRasterBand(2);
        in_bands[3] = in_sd_10m->GetRasterBand(3);
        in_bands[4] = in_sd_20m->GetRasterBand(1);
        in_bands[5] = in_sd_20m->GetRasterBand(2);
        in_bands[6] = in_sd_20m->GetRasterBand(3);
        in_bands[7] = in_sd_10m->GetRasterBand(4);
        in_bands[8] = in_sd_20m->GetRasterBand(4);
        in_bands[9] = in_sd_60m->GetRasterBand(2);
        in_bands[10] = in_sd_60m->GetRasterBand(3);
        in_bands[11] = in_sd_20m->GetRasterBand(5);
        in_bands[12] = in_sd_20m->GetRasterBand(6);




        for (uint16_t b=0; b<13; ++b) {
            if(out_vrt->AddBand(in_bands[b]->GetRasterDataType(), NULL) != CE_None) {
                std::cout << "ERROR: Cannot add band to VRT." << std::endl;
                return 1;
            }

            GDALRasterBand* x = out_vrt->GetRasterBand(out_vrt->GetRasterCount());

            char band_xml[10000];
            sprintf( band_xml,
                     "<SimpleSource resampling=\"%s\">"
                             "  <SourceFilename>%s</SourceFilename>"
                             "  <SourceBand>%d</SourceBand>"
                             "  <SourceProperties RasterXSize=\"%d\" RasterYSize=\"%d\" DataType=\"%s\" />"
                             "  <SrcRect xOff=\"0\" yOff=\"0\" xSize=\"%d\" ySize=\"%d\" />\n"
                             "  <DstRect xOff=\"0\" yOff=\"0\" xSize=\"%d\" ySize=\"%d\" />"
                             "  </SimpleSource>",
                     "nearest",
                     in_bands[b]->GetDataset()->GetDescription(),
                     in_bands[b]->GetBand(),
                     in_bands[b]->GetDataset()->GetRasterXSize(),
                     in_bands[b]->GetDataset()->GetRasterYSize(),
                     GDALGetDataTypeName(in_bands[b]->GetRasterDataType()),
                     in_bands[b]->GetDataset()->GetRasterXSize(),
                     in_bands[b]->GetDataset()->GetRasterYSize(),
                     in_sd_10m->GetRasterXSize(),
                     in_sd_10m->GetRasterXSize()
            );
            x->SetMetadataItem( "source_0", band_xml, "new_vrt_sources" );
        }


        GDALSetProjection(out_vrt, wkt_10m.c_str());
        double a[6];
        in_sd_10m->GetGeoTransform(a);
        GDALSetGeoTransform(out_vrt,a);




        if (!args.overview_levels.empty()) {

            if (args.verbose) {
                std::cout << "... building overviews (" << is + 1 << "/" << scene_sds.size() << ")" << std::endl;
            }
            // TODO generate overviews for in_sd_10m, in_sd_20m, and in_sd_60m
            std::vector<int> levels;
            for (uint16_t io=0; io<args.overview_levels.size(); ++io) {
                levels.push_back(args.overview_levels[io]);
            }

            if (in_sd_10m->BuildOverviews("NEAREST", args.overview_levels.size(), levels.data(), 0, NULL, GDALDummyProgress, NULL) != CE_None) {
                std::cout << "WARNING: generation of overviews failed." << std::endl; // TODO: make resampling variable
            }

            if (in_sd_20m->BuildOverviews("NEAREST", args.overview_levels.size(), levels.data(), 0, NULL, GDALDummyProgress, NULL) != CE_None) {
                std::cout << "WARNING: generation of overviews failed." << std::endl; // TODO: make resampling variable
            }

            if (in_sd_60m->BuildOverviews("NEAREST", args.overview_levels.size(), levels.data(), 0, NULL, GDALDummyProgress, NULL) != CE_None) {
                std::cout << "WARNING: generation of overviews failed." << std::endl; // TODO: make resampling variable
            }

        }






        GDALClose((GDALDatasetH)in_sd_10m);
        GDALClose((GDALDatasetH)in_sd_20m);
        GDALClose((GDALDatasetH)in_sd_60m);





        GDALClose((GDALDatasetH)out_vrt);

        ++it;
        ++is;
    }


}






int main(int ac, char *av[]) {






    // input args: input file, output directory, output format,  zoom levels, resampling,
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;


    po::options_description desc("Input arguments");
    desc.add_options()
            ("help,h", "print usage")
            ("overview,o", po::value<std::string>()->default_value(""), "overview levels which will be generated for input images")
            ("cog", "convert input Sentinel 2 subdatasets to cloud-optimized GeoTIFFs ")
            ("verbose,v", "verbose output")
            ("input-file", po::value<std::string>())
            ("output-file", po::value<std::string>());



    po::positional_options_description p;
    p.add("input-file", 1);
    p.add("output-file", 1);
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

    args.cog = false;
    if (vm.count("cog")) {
        args.cog = true;
    }




    if (!vm.count("input-file")) {
        std::cout << "ERROR: input file missing" << std::endl;
        return 1;
    }
    args.input_file = vm["input-file"].as<std::string>();


    if (!vm.count("output-file")) {
        std::cout << "ERROR: output file missing" << std::endl;
        return 1;
    }
    args.output_file = vm["output-file"].as<std::string>();


    args.overview_levels.clear();
    if (vm.count("overview") > 0) {
        std::istringstream str(vm["overview"].as<std::string>());
        std::vector<std::string> tokens{std::istream_iterator<std::string>{str}, std::istream_iterator<std::string>{}};
        std::vector<uint16_t> vec;
        for (uint8_t i = 0; i < tokens.size(); ++i) {
            vec.push_back((uint16_t) std::stoi(tokens[i]));
        }
        args.overview_levels = vec;
    }






    GDALAllRegister();

    create_vrt();




    return 0;
}
