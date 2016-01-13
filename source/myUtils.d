module myUtils;
import std.stdio;
import std.string;

void outputCSV3D(T)(T value, string filename, string spliter = ",", bool newfile = true) {
    if(newfile) {
        auto output = File(filename, "w");

        for(int i = 0; i < value.length; i++) {
            for(int j = 0; j < value.length; j++) {
                output.writeln(i, spliter, j, spliter, value[i][j]);
            }
            output.writeln();
        }
    } else {
        auto output = File(filename, "a");
        output.writeln("\n\n## hogehoge");

        foreach(int i, elems; value) {
            foreach(uint j, elem; elems) { 
                output.writeln(i, spliter, j, spliter, value[i][j]);
            }
            output.writeln();
        }
    }
}

void outputCSV2D(T)(T value, string filename, string spliter = ",", bool newfile = true) {
    if(newfile) {
        auto output = File(filename, "w");

        for(int i = 0; i < value.length; i++) {
            for(int j = 0; j < value.length; j++) {
                output.writeln(i, spliter, j, spliter, value[i][j]);
            }
            output.writeln();
        }
    } else {
        auto output = File(filename, "a");
        output.writeln("\n\n## hogehoge");

        foreach(int i, elems; value) {
            foreach(uint j, elem; elems) { 
                output.writeln(i, spliter, j, spliter, value[i][j]);
            }
            output.writeln();
        }
    }
}

void outputCSV1d(T)(T value, string filename, string spliter = ",", bool newfile = true) {
    if(newfile) {
        auto output = File(filename, "w");

        for(int i = 0; i < value.length; i++) {
            output.writeln(i, spliter, value[i]);
        }
    } else {
        auto output = File(filename, "a");
        output.writeln("\n\n## hogehoge");

        for(int i = 0; i < value.length; i++) {
            output.writeln(i, spliter, value[i]);
        }
    }
}

void print3D(T)(T data, const int nx, const int ny, const int nz){
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            for(int k = 0; k < nz; k++){
                write(data[i][j][k]);
                write(" ");
            }
            writeln();
        }
        writefln("----------", i);
    }

}
