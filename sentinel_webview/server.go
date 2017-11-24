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



package main

import (
	"net/http"
	"log"
	"github.com/gorilla/mux"
	"time"
	"os/exec"
	"strings"
	"encoding/json"
	"io/ioutil"
)



const EO_WEBTILE_PATH = "eo_webtile"
const TILE_DIR= "/tmp/tiles"  // directory, where PNG tiles will be stored
const SCENE_DIR= "data/" // directory, where S2 VRTs are stored
const CATALOG_FILE = "www/catalog.json"



type Scene struct {
	id    string
	path     string
	date	string
	footprint string
}




func main() {
	r := mux.NewRouter()
	//r.HandleFunc("/", req_root).Methods("GET")
	r.PathPrefix("/www/").Handler(http.StripPrefix("/www/", http.FileServer(http.Dir("www"))))
	r.HandleFunc("/tms/{scene_id}/{z}/{x}/{y}", serve_tile)
	r.HandleFunc("/search", search_data)
	srv := &http.Server{
		Handler:      r,
		Addr:         "localhost:8080",
		WriteTimeout: 10 * time.Second,
		ReadTimeout:  10 * time.Second,
	}
	log.Fatal(srv.ListenAndServe())
}




func serve_tile(w http.ResponseWriter, r *http.Request) {


	path_vars := mux.Vars(r)
	var bands string
	bands = "1 2 3" // default is first three bands

	query_params := r.URL.Query()
	if (len(query_params["bands"]) > 0) {
		s := strings.Split(query_params["bands"][0], ",")
		bands = s[0]
		for i := 1; i<len(s); i++ {
			bands += " " + s[i]
		}
	}

	var min string
	min = ""
	if (len(query_params["min"]) > 0) {
		s := strings.Split(query_params["min"][0], ",")
		min = s[0]
		for i := 1; i<len(s); i++ {
			min += " " + s[i]
		}
	}

	var max string
	max = ""
	if (len(query_params["max"]) > 0) {
		s := strings.Split(query_params["max"][0], ",")
		max = s[0]
		for i := 1; i<len(s); i++ {
			max += " " + s[i]
		}
	}

	var args []string

	args = append(args, "-b", bands)
	args = append(args, "-z", path_vars["z"])
	args = append(args, "-x", path_vars["x"])
	args = append(args, "-y", path_vars["y"])
	if (len(min) > 0) {
		args = append(args, "--min", min)
	}
	if (len(max) > 0) {
		args = append(args, "--max", max)
	}


	args = append(args, SCENE_DIR + "/" + path_vars["scene_id"] + ".vrt")
	args = append(args,  TILE_DIR)

	// TODO: check input query parameters and path variables for validity
	cmd := exec.Command(EO_WEBTILE_PATH, args...)
	err := cmd.Run()
	if err != nil {
		log.Fatal(err)
	}
	http.ServeFile(w, r, TILE_DIR + "/" + path_vars["z"] + "_" + path_vars["x"] + "_" + path_vars["y"] + ".png")
}




func search_data(w http.ResponseWriter, r *http.Request) {

	dat, err := ioutil.ReadFile(CATALOG_FILE)
	if err != nil {
		log.Fatal(err)
	}

	slist := make([]Scene,0)
	err = json.Unmarshal(dat, &slist)
	if err != nil {
		log.Fatal(err)
	}

	w.Header().Set("Content-Type", "application/json")
	w.Write(dat)
}





