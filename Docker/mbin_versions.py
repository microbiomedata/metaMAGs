#! /usr/bin/env python3
#njvarghese

def run():
    versions = {
        "mbin.py" : "0.5",
        "metabat2" : "2.15",
        "checkm-genome" : "1.2.1",
        "gtdb-tk" : "2.1.1",
        "hmmer" : "3.3.2",
        "prodigal" : "2.6.3",
        "pplacer" : "1.1.alpha19",
        "fasttree" : "2.1.11",
        "fastANI" : "1.33",
        "mash" : "2.3",
        "sqlite" : "3.39.2",
        "Python" : "3.9.12"
    }

    for name in versions.keys():
        print(name + ' : ' + versions[name])

if __name__ == '__main__':
    run()