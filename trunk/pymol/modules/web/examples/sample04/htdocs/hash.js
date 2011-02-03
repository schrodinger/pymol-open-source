function Hash(){
    for( var i=0; i < arguments.length; i++ )
        for( n in arguments[i] )
            if( arguments[i].hasOwnProperty(n) )
                this[n] = arguments[i][n];
}

    // Hash.version = 1.00;    // Original version
    // Hash.version = 1.01;    // Added ability to initialize in the constructor
    // Hash.version = 1.02;    // Fixed document bug that showed a non-working example (thanks mareks)
    //Hash.version = 1.03;    // Removed returning this from the constructor (thanks em-dash)
    Hash.version = 1.04;    // Missed some 'var' declarations (thanks Twey)


    Hash.prototype = new Object();

    Hash.prototype.keys = function(){
        var rv = [];
        for( var n in this )
            if( this.hasOwnProperty(n) )
                rv.push(n);
        return rv;
    }

    Hash.prototype.length = function(){
        return this.keys().length();
    }

    Hash.prototype.values = function(){
        var rv = [];
        for( var n in this )
            if( this.hasOwnProperty(n) )
                rv.push(this[n]);
        return rv;
    }

    Hash.prototype.slice = function(){
        var rv = [];
        for( var i = 0; i < arguments.length; i++ )
            rv.push(
                ( this.hasOwnProperty( arguments[i] ) )
                    ? this[arguments[i]]
                    : undefined
            );
        return rv;
    }

    Hash.prototype.concat = function(){
        for( var i = 0; i < arguments.length; i++ )
            for( var n in arguments[i] )
                if( arguments[i].hasOwnProperty(n) )
                    this[n] = arguments[i][n];
        return this;
    }
