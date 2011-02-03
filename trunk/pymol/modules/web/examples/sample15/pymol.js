// Copyright (C) Schrodinger, LLC, New York, NY. 
// All rights reserved.  

var pymol_fn_counter = 0;

function PyMOL(host, port, bufferMode, prefix) {

    this.Path = "/apply?_json="; // URL path to pass json args to pymol

    // now using _underscore_names for internal attributes & methods

    this._cmd_buffer = [];

    // do we throw up visible alerts when exceptions occur?

    this._alerts = true;

    this._host = host;
    this._port = port;

    this._parseBufferMode = function (bufferMode) {
        if ((bufferMode != undefined) && (bufferMode == 'on')) {
            return 'on';
        } else {
            return 'off';
        }
    }
    this._bufferMode = this._parseBufferMode(bufferMode);

    if(prefix == undefined) {
        this._prefix = 'pymol'; // default remote object name
    } else {
        this._prefix = prefix;
    }

    try {
        xmlhttp = new XMLHttpRequest();
    } catch(e) {
        xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }

    this.setBufferMode = function(bufferMode) {
        this._bufferMode = this._parseBufferMode(bufferMode);
        if (this._bufferMode == 'off') {
            this.flush();
        }
    }

    this.flush = function(callback) {
        if (this._cmd_buffer.length > 0) {
            var result = this._json('[' + this._cmd_buffer.join(',') + ']', callback);
            this._cmd_buffer.length = 0;
            return result;
        }
    }
    this.getBufferJSON = function() {
        return ('[' + this._cmd_buffer.join(',') + ']');
    }
    this.getBufferURL = function() {
        return (this.Path + '[' + this._cmd_buffer.join(',') + ']');
    }

    this._send_ajax = function(pypath, callback) {
        if (host == null) {
            myurl = pypath;
        } else {
            myurl = "http://" + host + ":" + port + pypath;
        }
        if(myurl.length>2000) { /* some broswers can't handle long URLs */
            var part = myurl.split("?",2);
            short_url = part[0];
            long_param = part[1];
            if (callback) {
                xmlhttp.open("POST", short_url, true);
                xmlhttp.onreadystatechange = callback;
            } else {
                xmlhttp.open("POST", short_url, false);
            }
            xmlhttp.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
            xmlhttp.setRequestHeader("Content-length", long_param.length);
            xmlhttp.setRequestHeader('Accept', 'text/json');
            xmlhttp.send(long_param);
        } else {
            if (callback) {
                xmlhttp.open("GET", myurl, true);
                xmlhttp.onreadystatechange = callback;
            } else {
                xmlhttp.open("GET", myurl, false);
            }
            xmlhttp.setRequestHeader('Accept', 'text/json');
            xmlhttp.send(null);
        }
        if (callback) {
        } else {
            if (xmlhttp.status == 500) {
                alert("PyMOL Exception:\n" + eval('(' + xmlhttp.responseText + ')').join("\n"));
                return null;
            } 
            return eval('(' + xmlhttp.responseText + ')');
        }
        return false;
    }
    
    this._send_cross_script = function(pypath, callback) {
        if (host == null) {
            myurl = pypath;
        } else {
            myurl = "http://" + host + ":" + port + pypath;
        }
        if (callback == undefined) {
            //myurl +=  "&_callback=alert"
            myurl +=  "&_callback=void"
        } else {
            // IE does not have callback.name, so we invent a unique one
            // and store it in the "global" namespace, window.  In this way,
            // the javascript callback from our cross-domain script hack
            // contains a javascript function name known to this page/window.
            if (callback.name == undefined) {
               callback.name = '_fn' + pymol_fn_counter++;
               window[ callback.name ] = callback;
            } 
            myurl +=  "&_callback=" + callback.name;
        }
        var head = document.getElementsByTagName("head")[0];
        var script = document.createElement("script");
        script.src = myurl;
        head.appendChild(script);
        return true;
    }

    this._send = function(pypath, callback) {
        if ( (host == null) || 
             ((host == document.domain) && (port == document.location.port)) ) {
            return this._send_ajax(pypath, callback);
        } else {
            return this._send_cross_script(pypath, callback);
        }
    }
    
    this._handle_response = function(status, text) {
        response = JSON.parse(text);
        if(status == 200) { // normal result
            return response;
        } else { // some kind of error condition
            if(this._alerts) {
                alert(response.join("\n"));
            }
            return;
        }
    }

    this._result = function(e) {
        if ((host == document.domain) || (host == null)) {
            if (typeof e == "object") {
                // asynchronous ajax was used
                if (xmlhttp.readyState == 4) {
                    return this._handle_response(xmlhttp.status, xmlhttp.responseText);
                    return;
                }
            }
            // "synchronous" (no-callback) ajax was used
            if(xmlhttp.responseText) {
                return this._handle_response(xmlhttp.status, xmlhttp.responseText);
            } else {
                return;
            }
        } 
        // e was provided by cross-script callback, so just return it
        return e;
    }

    // for private use
    this._json = function(jcmd,callback) {
        return this._send(this.Path+jcmd,callback);
    }
    // for public use;  args switched to match style of pymol.cmd calls
    this.sendJSON = function(anyargs) {
        if (typeof arguments[0] == 'function') {
            //this.sendJSON = function(callback, jcmd) {
            callback = arguments[0];
            jcmd     = arguments[1];
        } else {
            //this.sendJSON = function(jcmd) {
            callback = null;
            jcmd     = arguments[0];
        }
        return this._send(this.Path+jcmd,callback);
    }

    this._apply = function(name, args, kwds, callback) {
        // note: javascript's 'this' can refer to either pymol or cmd
        if (name.substring(0,1) == '.') { // .cmd.method -> pymol.cmd.method
            name = this._prefix + name;
        }
        var mypath = this.Path;
        var myargs = '["' + name + '"';
        if (args != null)  {
            myargs = myargs + "," + JSON.stringify(args);
            if(kwds != null) {
                myargs = myargs + "," + JSON.stringify(kwds);
            }
        } else if(kwds != null) {
            myargs = myargs + ",[]," + JSON.stringify(kwds)
        }
        myargs = myargs + "]";

        //        document.getElementById('debug').innerHTML = mypath;

        if (this._bufferMode == 'on') {
            this._cmd_buffer.push(myargs);
            if ( callback ) {
                this.flush(callback);
            } else {
                return;
            }
        } else {
            return this._send(mypath+myargs,callback);
        }
    }

    this.apply_cmd = function(name, args, kwds, callback) {
        this._apply(".cmd." + name, args, kwds, callback);
    }

    this.cmd = { // cmd is a public attribute of PyMOL instances
        
        _dispatch: function(anyargs) {
           //  name is always first argument
           var name = ".cmd." + arguments[0];
           // callback function is optional second argument
           var callback = undefined;
           var argstart = 1;
           if (typeof arguments[1] == 'function') {
             callback = arguments[1];
             argstart = 2;
           }
           return this._outer._apply(name,Array.prototype.slice.apply(arguments).slice(argstart),null,callback);
        },

        // BEGIN MACHINE-GENERATED CODE
        // python webapi.py

        load: function() {
            return this._dispatch.apply(this,
                ["load"].concat(Array.prototype.slice.apply(arguments)));
        },
        load_traj: function() {
            return this._dispatch.apply(this,
                ["load_traj"].concat(Array.prototype.slice.apply(arguments)));
        },
        load_png: function() {
            return this._dispatch.apply(this,
                ["load_png"].concat(Array.prototype.slice.apply(arguments)));
        },
        fragment: function() {
            return this._dispatch.apply(this,
                ["fragment"].concat(Array.prototype.slice.apply(arguments)));
        },
        fetch: function() {
            return this._dispatch.apply(this,
                ["fetch"].concat(Array.prototype.slice.apply(arguments)));
        },
        read_mmodstr: function() {
            return this._dispatch.apply(this,
                ["read_mmodstr"].concat(Array.prototype.slice.apply(arguments)));
        },
        read_molstr: function() {
            return this._dispatch.apply(this,
                ["read_molstr"].concat(Array.prototype.slice.apply(arguments)));
        },
        read_sdfstr: function() {
            return this._dispatch.apply(this,
                ["read_sdfstr"].concat(Array.prototype.slice.apply(arguments)));
        },
        read_pdbstr: function() {
            return this._dispatch.apply(this,
                ["read_pdbstr"].concat(Array.prototype.slice.apply(arguments)));
        },
        read_xplorstr: function() {
            return this._dispatch.apply(this,
                ["read_xplorstr"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_pdbstr: function() {
            return this._dispatch.apply(this,
                ["get_pdbstr"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_fastastr: function() {
            return this._dispatch.apply(this,
                ["get_fastastr"].concat(Array.prototype.slice.apply(arguments)));
        },
        copy: function() {
            return this._dispatch.apply(this,
                ["copy"].concat(Array.prototype.slice.apply(arguments)));
        },
        create: function() {
            return this._dispatch.apply(this,
                ["create"].concat(Array.prototype.slice.apply(arguments)));
        },
        extract: function() {
            return this._dispatch.apply(this,
                ["extract"].concat(Array.prototype.slice.apply(arguments)));
        },
        split_states: function() {
            return this._dispatch.apply(this,
                ["split_states"].concat(Array.prototype.slice.apply(arguments)));
        },
        symexp: function() {
            return this._dispatch.apply(this,
                ["symexp"].concat(Array.prototype.slice.apply(arguments)));
        },
        ramp_new: function() {
            return this._dispatch.apply(this,
                ["ramp_new"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_name: function() {
            return this._dispatch.apply(this,
                ["set_name"].concat(Array.prototype.slice.apply(arguments)));
        },
        map_new: function() {
            return this._dispatch.apply(this,
                ["map_new"].concat(Array.prototype.slice.apply(arguments)));
        },
        map_set: function() {
            return this._dispatch.apply(this,
                ["map_set"].concat(Array.prototype.slice.apply(arguments)));
        },
        map_set_border: function() {
            return this._dispatch.apply(this,
                ["map_set_border"].concat(Array.prototype.slice.apply(arguments)));
        },
        map_double: function() {
            return this._dispatch.apply(this,
                ["map_double"].concat(Array.prototype.slice.apply(arguments)));
        },
        map_halve: function() {
            return this._dispatch.apply(this,
                ["map_halve"].concat(Array.prototype.slice.apply(arguments)));
        },
        map_trim: function() {
            return this._dispatch.apply(this,
                ["map_trim"].concat(Array.prototype.slice.apply(arguments)));
        },
        isodot: function() {
            return this._dispatch.apply(this,
                ["isodot"].concat(Array.prototype.slice.apply(arguments)));
        },
        isolevel: function() {
            return this._dispatch.apply(this,
                ["isolevel"].concat(Array.prototype.slice.apply(arguments)));
        },
        isomesh: function() {
            return this._dispatch.apply(this,
                ["isomesh"].concat(Array.prototype.slice.apply(arguments)));
        },
        isosurface: function() {
            return this._dispatch.apply(this,
                ["isosurface"].concat(Array.prototype.slice.apply(arguments)));
        },
        slice_new: function() {
            return this._dispatch.apply(this,
                ["slice_new"].concat(Array.prototype.slice.apply(arguments)));
        },
        gradient: function() {
            return this._dispatch.apply(this,
                ["gradient"].concat(Array.prototype.slice.apply(arguments)));
        },
        ungroup: function() {
            return this._dispatch.apply(this,
                ["ungroup"].concat(Array.prototype.slice.apply(arguments)));
        },
        group: function() {
            return this._dispatch.apply(this,
                ["group"].concat(Array.prototype.slice.apply(arguments)));
        },
        pseudoatom: function() {
            return this._dispatch.apply(this,
                ["pseudoatom"].concat(Array.prototype.slice.apply(arguments)));
        },
        fab: function() {
            return this._dispatch.apply(this,
                ["fab"].concat(Array.prototype.slice.apply(arguments)));
        },
        enable: function() {
            return this._dispatch.apply(this,
                ["enable"].concat(Array.prototype.slice.apply(arguments)));
        },
        disable: function() {
            return this._dispatch.apply(this,
                ["disable"].concat(Array.prototype.slice.apply(arguments)));
        },
        delete_: function() {
            return this._dispatch.apply(this,
                ["delete_"].concat(Array.prototype.slice.apply(arguments)));
        },
        reinitialize: function() {
            return this._dispatch.apply(this,
                ["reinitialize"].concat(Array.prototype.slice.apply(arguments)));
        },
        deselect: function() {
            return this._dispatch.apply(this,
                ["deselect"].concat(Array.prototype.slice.apply(arguments)));
        },
        select: function() {
            return this._dispatch.apply(this,
                ["select"].concat(Array.prototype.slice.apply(arguments)));
        },
        indicate: function() {
            return this._dispatch.apply(this,
                ["indicate"].concat(Array.prototype.slice.apply(arguments)));
        },
        select_list: function() {
            return this._dispatch.apply(this,
                ["select_list"].concat(Array.prototype.slice.apply(arguments)));
        },
        pop: function() {
            return this._dispatch.apply(this,
                ["pop"].concat(Array.prototype.slice.apply(arguments)));
        },
        angle: function() {
            return this._dispatch.apply(this,
                ["angle"].concat(Array.prototype.slice.apply(arguments)));
        },
        dihedral: function() {
            return this._dispatch.apply(this,
                ["dihedral"].concat(Array.prototype.slice.apply(arguments)));
        },
        dist: function() {
            return this._dispatch.apply(this,
                ["dist"].concat(Array.prototype.slice.apply(arguments)));
        },
        distance: function() {
            return this._dispatch.apply(this,
                ["distance"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_angle: function() {
            return this._dispatch.apply(this,
                ["get_angle"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_dihedral: function() {
            return this._dispatch.apply(this,
                ["get_dihedral"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_distance: function() {
            return this._dispatch.apply(this,
                ["get_distance"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_area: function() {
            return this._dispatch.apply(this,
                ["get_area"].concat(Array.prototype.slice.apply(arguments)));
        },
        color: function() {
            return this._dispatch.apply(this,
                ["color"].concat(Array.prototype.slice.apply(arguments)));
        },
        bg_color: function() {
            return this._dispatch.apply(this,
                ["bg_color"].concat(Array.prototype.slice.apply(arguments)));
        },
        rebuild: function() {
            return this._dispatch.apply(this,
                ["rebuild"].concat(Array.prototype.slice.apply(arguments)));
        },
        refresh: function() {
            return this._dispatch.apply(this,
                ["refresh"].concat(Array.prototype.slice.apply(arguments)));
        },
        recolor: function() {
            return this._dispatch.apply(this,
                ["recolor"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_color: function() {
            return this._dispatch.apply(this,
                ["set_color"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_object_color: function() {
            return this._dispatch.apply(this,
                ["set_object_color"].concat(Array.prototype.slice.apply(arguments)));
        },
        show: function() {
            return this._dispatch.apply(this,
                ["show"].concat(Array.prototype.slice.apply(arguments)));
        },
        show_as: function() {
            return this._dispatch.apply(this,
                ["show_as"].concat(Array.prototype.slice.apply(arguments)));
        },
        hide: function() {
            return this._dispatch.apply(this,
                ["hide"].concat(Array.prototype.slice.apply(arguments)));
        },
        cartoon: function() {
            return this._dispatch.apply(this,
                ["cartoon"].concat(Array.prototype.slice.apply(arguments)));
        },
        spectrum: function() {
            return this._dispatch.apply(this,
                ["spectrum"].concat(Array.prototype.slice.apply(arguments)));
        },
        center: function() {
            return this._dispatch.apply(this,
                ["center"].concat(Array.prototype.slice.apply(arguments)));
        },
        zoom: function() {
            return this._dispatch.apply(this,
                ["zoom"].concat(Array.prototype.slice.apply(arguments)));
        },
        reset: function() {
            return this._dispatch.apply(this,
                ["reset"].concat(Array.prototype.slice.apply(arguments)));
        },
        clip: function() {
            return this._dispatch.apply(this,
                ["clip"].concat(Array.prototype.slice.apply(arguments)));
        },
        orient: function() {
            return this._dispatch.apply(this,
                ["orient"].concat(Array.prototype.slice.apply(arguments)));
        },
        origin: function() {
            return this._dispatch.apply(this,
                ["origin"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_view: function() {
            return this._dispatch.apply(this,
                ["set_view"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_view: function() {
            return this._dispatch.apply(this,
                ["get_view"].concat(Array.prototype.slice.apply(arguments)));
        },
        move: function() {
            return this._dispatch.apply(this,
                ["move"].concat(Array.prototype.slice.apply(arguments)));
        },
        turn: function() {
            return this._dispatch.apply(this,
                ["turn"].concat(Array.prototype.slice.apply(arguments)));
        },
        rock: function() {
            return this._dispatch.apply(this,
                ["rock"].concat(Array.prototype.slice.apply(arguments)));
        },
        stereo: function() {
            return this._dispatch.apply(this,
                ["stereo"].concat(Array.prototype.slice.apply(arguments)));
        },
        get: function() {
            return this._dispatch.apply(this,
                ["get"].concat(Array.prototype.slice.apply(arguments)));
        },
        set: function() {
            return this._dispatch.apply(this,
                ["set"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_bond: function() {
            return this._dispatch.apply(this,
                ["set_bond"].concat(Array.prototype.slice.apply(arguments)));
        },
        unset: function() {
            return this._dispatch.apply(this,
                ["unset"].concat(Array.prototype.slice.apply(arguments)));
        },
        unset_bond: function() {
            return this._dispatch.apply(this,
                ["unset_bond"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_setting_boolean: function() {
            return this._dispatch.apply(this,
                ["get_setting_boolean"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_setting_int: function() {
            return this._dispatch.apply(this,
                ["get_setting_int"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_setting_float: function() {
            return this._dispatch.apply(this,
                ["get_setting_float"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_setting_legacy: function() {
            return this._dispatch.apply(this,
                ["get_setting_legacy"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_setting_tuple: function() {
            return this._dispatch.apply(this,
                ["get_setting_tuple"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_setting_text: function() {
            return this._dispatch.apply(this,
                ["get_setting_text"].concat(Array.prototype.slice.apply(arguments)));
        },
        window: function() {
            return this._dispatch.apply(this,
                ["window"].concat(Array.prototype.slice.apply(arguments)));
        },
        viewport: function() {
            return this._dispatch.apply(this,
                ["viewport"].concat(Array.prototype.slice.apply(arguments)));
        },
        full_screen: function() {
            return this._dispatch.apply(this,
                ["full_screen"].concat(Array.prototype.slice.apply(arguments)));
        },
        quit: function() {
            return this._dispatch.apply(this,
                ["quit"].concat(Array.prototype.slice.apply(arguments)));
        },
        draw: function() {
            return this._dispatch.apply(this,
                ["draw"].concat(Array.prototype.slice.apply(arguments)));
        },
        ray: function() {
            return this._dispatch.apply(this,
                ["ray"].concat(Array.prototype.slice.apply(arguments)));
        },
        align: function() {
            return this._dispatch.apply(this,
                ["align"].concat(Array.prototype.slice.apply(arguments)));
        },
        super_: function() {
            return this._dispatch.apply(this,
                ["super_"].concat(Array.prototype.slice.apply(arguments)));
        },
        fit: function() {
            return this._dispatch.apply(this,
                ["fit"].concat(Array.prototype.slice.apply(arguments)));
        },
        rms: function() {
            return this._dispatch.apply(this,
                ["rms"].concat(Array.prototype.slice.apply(arguments)));
        },
        rms_cur: function() {
            return this._dispatch.apply(this,
                ["rms_cur"].concat(Array.prototype.slice.apply(arguments)));
        },
        intra_fit: function() {
            return this._dispatch.apply(this,
                ["intra_fit"].concat(Array.prototype.slice.apply(arguments)));
        },
        intra_rms: function() {
            return this._dispatch.apply(this,
                ["intra_rms"].concat(Array.prototype.slice.apply(arguments)));
        },
        intra_rms_cur: function() {
            return this._dispatch.apply(this,
                ["intra_rms_cur"].concat(Array.prototype.slice.apply(arguments)));
        },
        pair_fit: function() {
            return this._dispatch.apply(this,
                ["pair_fit"].concat(Array.prototype.slice.apply(arguments)));
        },
        space: function() {
            return this._dispatch.apply(this,
                ["space"].concat(Array.prototype.slice.apply(arguments)));
        },
        order: function() {
            return this._dispatch.apply(this,
                ["order"].concat(Array.prototype.slice.apply(arguments)));
        },
        edit_mode: function() {
            return this._dispatch.apply(this,
                ["edit_mode"].concat(Array.prototype.slice.apply(arguments)));
        },
        button: function() {
            return this._dispatch.apply(this,
                ["button"].concat(Array.prototype.slice.apply(arguments)));
        },
        config_mouse: function() {
            return this._dispatch.apply(this,
                ["config_mouse"].concat(Array.prototype.slice.apply(arguments)));
        },
        mouse: function() {
            return this._dispatch.apply(this,
                ["mouse"].concat(Array.prototype.slice.apply(arguments)));
        },
        mask: function() {
            return this._dispatch.apply(this,
                ["mask"].concat(Array.prototype.slice.apply(arguments)));
        },
        unmask: function() {
            return this._dispatch.apply(this,
                ["unmask"].concat(Array.prototype.slice.apply(arguments)));
        },
        count_atoms: function() {
            return this._dispatch.apply(this,
                ["count_atoms"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_chains: function() {
            return this._dispatch.apply(this,
                ["get_chains"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_color_index: function() {
            return this._dispatch.apply(this,
                ["get_color_index"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_color_indices: function() {
            return this._dispatch.apply(this,
                ["get_color_indices"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_object_color_index: function() {
            return this._dispatch.apply(this,
                ["get_object_color_index"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_object_list: function() {
            return this._dispatch.apply(this,
                ["get_object_list"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_color_tuple: function() {
            return this._dispatch.apply(this,
                ["get_color_tuple"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_atom_coords: function() {
            return this._dispatch.apply(this,
                ["get_atom_coords"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_extent: function() {
            return this._dispatch.apply(this,
                ["get_extent"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_names: function() {
            return this._dispatch.apply(this,
                ["get_names"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_names_of_type: function() {
            return this._dispatch.apply(this,
                ["get_names_of_type"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_legal_name: function() {
            return this._dispatch.apply(this,
                ["get_legal_name"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_unused_name: function() {
            return this._dispatch.apply(this,
                ["get_unused_name"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_object_matrix: function() {
            return this._dispatch.apply(this,
                ["get_object_matrix"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_phipsi: function() {
            return this._dispatch.apply(this,
                ["get_phipsi"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_position: function() {
            return this._dispatch.apply(this,
                ["get_position"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_raw_alignment: function() {
            return this._dispatch.apply(this,
                ["get_raw_alignment"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_renderer: function() {
            return this._dispatch.apply(this,
                ["get_renderer"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_symmetry: function() {
            return this._dispatch.apply(this,
                ["get_symmetry"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_title: function() {
            return this._dispatch.apply(this,
                ["get_title"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_type: function() {
            return this._dispatch.apply(this,
                ["get_type"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_version: function() {
            return this._dispatch.apply(this,
                ["get_version"].concat(Array.prototype.slice.apply(arguments)));
        },
        id_atom: function() {
            return this._dispatch.apply(this,
                ["id_atom"].concat(Array.prototype.slice.apply(arguments)));
        },
        identify: function() {
            return this._dispatch.apply(this,
                ["identify"].concat(Array.prototype.slice.apply(arguments)));
        },
        index: function() {
            return this._dispatch.apply(this,
                ["index"].concat(Array.prototype.slice.apply(arguments)));
        },
        phi_psi: function() {
            return this._dispatch.apply(this,
                ["phi_psi"].concat(Array.prototype.slice.apply(arguments)));
        },
        matrix_copy: function() {
            return this._dispatch.apply(this,
                ["matrix_copy"].concat(Array.prototype.slice.apply(arguments)));
        },
        matrix_reset: function() {
            return this._dispatch.apply(this,
                ["matrix_reset"].concat(Array.prototype.slice.apply(arguments)));
        },
        rotate: function() {
            return this._dispatch.apply(this,
                ["rotate"].concat(Array.prototype.slice.apply(arguments)));
        },
        translate: function() {
            return this._dispatch.apply(this,
                ["translate"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_object_ttt: function() {
            return this._dispatch.apply(this,
                ["set_object_ttt"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_dihedral: function() {
            return this._dispatch.apply(this,
                ["set_dihedral"].concat(Array.prototype.slice.apply(arguments)));
        },
        transform_object: function() {
            return this._dispatch.apply(this,
                ["transform_object"].concat(Array.prototype.slice.apply(arguments)));
        },
        transform_selection: function() {
            return this._dispatch.apply(this,
                ["transform_selection"].concat(Array.prototype.slice.apply(arguments)));
        },
        translate_atom: function() {
            return this._dispatch.apply(this,
                ["translate_atom"].concat(Array.prototype.slice.apply(arguments)));
        },
        update: function() {
            return this._dispatch.apply(this,
                ["update"].concat(Array.prototype.slice.apply(arguments)));
        },
        attach: function() {
            return this._dispatch.apply(this,
                ["attach"].concat(Array.prototype.slice.apply(arguments)));
        },
        bond: function() {
            return this._dispatch.apply(this,
                ["bond"].concat(Array.prototype.slice.apply(arguments)));
        },
        unbond: function() {
            return this._dispatch.apply(this,
                ["unbond"].concat(Array.prototype.slice.apply(arguments)));
        },
        cycle_valence: function() {
            return this._dispatch.apply(this,
                ["cycle_valence"].concat(Array.prototype.slice.apply(arguments)));
        },
        drag: function() {
            return this._dispatch.apply(this,
                ["drag"].concat(Array.prototype.slice.apply(arguments)));
        },
        dss: function() {
            return this._dispatch.apply(this,
                ["dss"].concat(Array.prototype.slice.apply(arguments)));
        },
        edit: function() {
            return this._dispatch.apply(this,
                ["edit"].concat(Array.prototype.slice.apply(arguments)));
        },
        unpick: function() {
            return this._dispatch.apply(this,
                ["unpick"].concat(Array.prototype.slice.apply(arguments)));
        },
        fix_chemistry: function() {
            return this._dispatch.apply(this,
                ["fix_chemistry"].concat(Array.prototype.slice.apply(arguments)));
        },
        flag: function() {
            return this._dispatch.apply(this,
                ["flag"].concat(Array.prototype.slice.apply(arguments)));
        },
        fuse: function() {
            return this._dispatch.apply(this,
                ["fuse"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_editor_scheme: function() {
            return this._dispatch.apply(this,
                ["get_editor_scheme"].concat(Array.prototype.slice.apply(arguments)));
        },
        h_add: function() {
            return this._dispatch.apply(this,
                ["h_add"].concat(Array.prototype.slice.apply(arguments)));
        },
        h_fill: function() {
            return this._dispatch.apply(this,
                ["h_fill"].concat(Array.prototype.slice.apply(arguments)));
        },
        h_fix: function() {
            return this._dispatch.apply(this,
                ["h_fix"].concat(Array.prototype.slice.apply(arguments)));
        },
        invert: function() {
            return this._dispatch.apply(this,
                ["invert"].concat(Array.prototype.slice.apply(arguments)));
        },
        torsion: function() {
            return this._dispatch.apply(this,
                ["torsion"].concat(Array.prototype.slice.apply(arguments)));
        },
        valence: function() {
            return this._dispatch.apply(this,
                ["valence"].concat(Array.prototype.slice.apply(arguments)));
        },
        clean: function() {
            return this._dispatch.apply(this,
                ["clean"].concat(Array.prototype.slice.apply(arguments)));
        },
        deprotect: function() {
            return this._dispatch.apply(this,
                ["deprotect"].concat(Array.prototype.slice.apply(arguments)));
        },
        protect: function() {
            return this._dispatch.apply(this,
                ["protect"].concat(Array.prototype.slice.apply(arguments)));
        },
        reference: function() {
            return this._dispatch.apply(this,
                ["reference"].concat(Array.prototype.slice.apply(arguments)));
        },
        remove: function() {
            return this._dispatch.apply(this,
                ["remove"].concat(Array.prototype.slice.apply(arguments)));
        },
        remove_picked: function() {
            return this._dispatch.apply(this,
                ["remove_picked"].concat(Array.prototype.slice.apply(arguments)));
        },
        rename: function() {
            return this._dispatch.apply(this,
                ["rename"].concat(Array.prototype.slice.apply(arguments)));
        },
        replace: function() {
            return this._dispatch.apply(this,
                ["replace"].concat(Array.prototype.slice.apply(arguments)));
        },
        sculpt_purge: function() {
            return this._dispatch.apply(this,
                ["sculpt_purge"].concat(Array.prototype.slice.apply(arguments)));
        },
        sculpt_deactivate: function() {
            return this._dispatch.apply(this,
                ["sculpt_deactivate"].concat(Array.prototype.slice.apply(arguments)));
        },
        sculpt_activate: function() {
            return this._dispatch.apply(this,
                ["sculpt_activate"].concat(Array.prototype.slice.apply(arguments)));
        },
        sculpt_iterate: function() {
            return this._dispatch.apply(this,
                ["sculpt_iterate"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_geometry: function() {
            return this._dispatch.apply(this,
                ["set_geometry"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_symmetry: function() {
            return this._dispatch.apply(this,
                ["set_symmetry"].concat(Array.prototype.slice.apply(arguments)));
        },
        set_title: function() {
            return this._dispatch.apply(this,
                ["set_title"].concat(Array.prototype.slice.apply(arguments)));
        },
        smooth: function() {
            return this._dispatch.apply(this,
                ["smooth"].concat(Array.prototype.slice.apply(arguments)));
        },
        sort: function() {
            return this._dispatch.apply(this,
                ["sort"].concat(Array.prototype.slice.apply(arguments)));
        },
        undo: function() {
            return this._dispatch.apply(this,
                ["undo"].concat(Array.prototype.slice.apply(arguments)));
        },
        push_undo: function() {
            return this._dispatch.apply(this,
                ["push_undo"].concat(Array.prototype.slice.apply(arguments)));
        },
        redo: function() {
            return this._dispatch.apply(this,
                ["redo"].concat(Array.prototype.slice.apply(arguments)));
        },
        wizard: function() {
            return this._dispatch.apply(this,
                ["wizard"].concat(Array.prototype.slice.apply(arguments)));
        },
        replace_wizard: function() {
            return this._dispatch.apply(this,
                ["replace_wizard"].concat(Array.prototype.slice.apply(arguments)));
        },
        count_frames: function() {
            return this._dispatch.apply(this,
                ["count_frames"].concat(Array.prototype.slice.apply(arguments)));
        },
        count_states: function() {
            return this._dispatch.apply(this,
                ["count_states"].concat(Array.prototype.slice.apply(arguments)));
        },
        mset: function() {
            return this._dispatch.apply(this,
                ["mset"].concat(Array.prototype.slice.apply(arguments)));
        },
        madd: function() {
            return this._dispatch.apply(this,
                ["madd"].concat(Array.prototype.slice.apply(arguments)));
        },
        mclear: function() {
            return this._dispatch.apply(this,
                ["mclear"].concat(Array.prototype.slice.apply(arguments)));
        },
        mmatrix: function() {
            return this._dispatch.apply(this,
                ["mmatrix"].concat(Array.prototype.slice.apply(arguments)));
        },
        mdump: function() {
            return this._dispatch.apply(this,
                ["mdump"].concat(Array.prototype.slice.apply(arguments)));
        },
        mview: function() {
            return this._dispatch.apply(this,
                ["mview"].concat(Array.prototype.slice.apply(arguments)));
        },
        forward: function() {
            return this._dispatch.apply(this,
                ["forward"].concat(Array.prototype.slice.apply(arguments)));
        },
        backward: function() {
            return this._dispatch.apply(this,
                ["backward"].concat(Array.prototype.slice.apply(arguments)));
        },
        rewind: function() {
            return this._dispatch.apply(this,
                ["rewind"].concat(Array.prototype.slice.apply(arguments)));
        },
        middle: function() {
            return this._dispatch.apply(this,
                ["middle"].concat(Array.prototype.slice.apply(arguments)));
        },
        ending: function() {
            return this._dispatch.apply(this,
                ["ending"].concat(Array.prototype.slice.apply(arguments)));
        },
        mplay: function() {
            return this._dispatch.apply(this,
                ["mplay"].concat(Array.prototype.slice.apply(arguments)));
        },
        mtoggle: function() {
            return this._dispatch.apply(this,
                ["mtoggle"].concat(Array.prototype.slice.apply(arguments)));
        },
        mstop: function() {
            return this._dispatch.apply(this,
                ["mstop"].concat(Array.prototype.slice.apply(arguments)));
        },
        frame: function() {
            return this._dispatch.apply(this,
                ["frame"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_movie_playing: function() {
            return this._dispatch.apply(this,
                ["get_movie_playing"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_state: function() {
            return this._dispatch.apply(this,
                ["get_state"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_frame: function() {
            return this._dispatch.apply(this,
                ["get_frame"].concat(Array.prototype.slice.apply(arguments)));
        },
        view: function() {
            return this._dispatch.apply(this,
                ["view"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_scene_dict: function() {
            return this._dispatch.apply(this,
                ["get_scene_dict"].concat(Array.prototype.slice.apply(arguments)));
        },
        get_scene_list: function() {
            return this._dispatch.apply(this,
                ["get_scene_list"].concat(Array.prototype.slice.apply(arguments)));
        },
        scene: function() {
            return this._dispatch.apply(this,
                ["scene"].concat(Array.prototype.slice.apply(arguments)));
        },
        scene_order: function() {
            return this._dispatch.apply(this,
                ["scene_order"].concat(Array.prototype.slice.apply(arguments)));
        }

        // END MACHINE-GENERATED CODE

    }
    
    this.getattr = {
      dispatch: function(attrname) {
        if (attrname) {
          mypath = "/getattr/pymol." + attrname;
          this._outer._send(mypath);
        }
      }
      //      viewing: function() {
      //        this._dispatch("viewing");
      //      }
    }

    // enable inner pseudo-inner-instance to find the outer instance

    this.cmd._outer                = this
    this.getattr._outer            = this
}

// utility functions not used above, but of use to some users, one hopes
function parse_file(f) {
  var re = /[\\\/]/;
  fields = f.split(re);
  filename = fields[fields.length-1];
  dot = filename.lastIndexOf('.');
  if (dot > -1) {
    name = filename.slice(0,dot);
    ext    = filename.slice(dot);
  } else {
    name = filename;
    ext = "";
  }
  // Javascript 1.7+ only in Mozilla?
  // return [name, ext];
  return new Array(name, ext);
}
function validate_file(n,f) {
  if (f.value) {
    // Javascript 1.7+ only in Mozilla?
    // [name, ext] = parse_file(f.value);
    parsed = parse_file(f.value);
    name = parsed[0];
    ext = parsed[1];
    if (name) {
      n.value = name;
      //add_checkbox(name, "objects");
      return true;
    } else {
      alert('cannot extract name from file');
      return false;
    }
  } else {
    alert('no file selected');
    return false;
  }
}
