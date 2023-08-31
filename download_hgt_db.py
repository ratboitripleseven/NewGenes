import os
import subprocess

ROOT_FOLDER = "./temp_data"


LINK_TO_FILES = [
    "https://usuaris.tinet.cat/debb/HGT/aquae.d/aquae.doc",
    "https://usuaris.tinet.cat/debb/HGT/cmur.d/cmur.doc",
    "https://usuaris.tinet.cat/debb/HGT/ctra.d/ctra.doc",
    "https://usuaris.tinet.cat/debb/HGT/cpneu.d/cpneu.doc",
    "https://usuaris.tinet.cat/debb/HGT/cpneu2.d/cpneu2.doc",
    "https://usuaris.tinet.cat/debb/HGT/cpneu3.d/cpneu3.doc",
    "https://usuaris.tinet.cat/debb/HGT/btheta.d/btheta.doc",
    "https://usuaris.tinet.cat/debb/HGT/ctep.d/ctep.doc",
    "https://usuaris.tinet.cat/debb/HGT/synecho.d/synecho.doc",
    "https://usuaris.tinet.cat/debb/HGT/nsp.d/nsp.doc",
    "https://usuaris.tinet.cat/debb/HGT/cglu.d/cglu.doc",
    "https://usuaris.tinet.cat/debb/HGT/mtub.d/mtub.doc",
    "https://usuaris.tinet.cat/debb/HGT/mtub2.d/mtub2.doc",
    "https://usuaris.tinet.cat/debb/HGT/mlep.d/mlep.doc",
    "https://usuaris.tinet.cat/debb/HGT/scoel.d/scoel.doc",
    "https://usuaris.tinet.cat/debb/HGT/bsub.d/bsub.doc",
    "https://usuaris.tinet.cat/debb/HGT/bhal.d/bhal.doc",
    "https://usuaris.tinet.cat/debb/HGT/linno.d/linno.doc",
    "https://usuaris.tinet.cat/debb/HGT/lmono.d/lmono.doc",
    "https://usuaris.tinet.cat/debb/HGT/sau2.d/sau2.doc",
    "https://usuaris.tinet.cat/debb/HGT/sau1.d/sau1.doc",
    "https://usuaris.tinet.cat/debb/HGT/sau3.d/sau3.doc",
    "https://usuaris.tinet.cat/debb/HGT/llac.d/llac.doc",
    "https://usuaris.tinet.cat/debb/HGT/spyo.d/spyo.doc",
    "https://usuaris.tinet.cat/debb/HGT/spyo2.d/spyo2.doc",
    "https://usuaris.tinet.cat/debb/HGT/spneu1.d/spneu1.doc",
    "https://usuaris.tinet.cat/debb/HGT/spneu2.d/spneu2.doc",
    "https://usuaris.tinet.cat/debb/HGT/caceto.d/caceto.doc",
    "https://usuaris.tinet.cat/debb/HGT/cperf.d/cperf.doc",
    "https://usuaris.tinet.cat/debb/HGT/tteng.d/tteng.doc",
    "https://usuaris.tinet.cat/debb/HGT/mgen.d/mgen.doc",
    "https://usuaris.tinet.cat/debb/HGT/mpneu.d/mpneu.doc",
    "https://usuaris.tinet.cat/debb/HGT/mpul.d/mpul.doc",
    "https://usuaris.tinet.cat/debb/HGT/uure.d/uure.doc",
    "https://usuaris.tinet.cat/debb/HGT/fnucl.d/fnucl.doc",
    "https://usuaris.tinet.cat/debb/HGT/ccres.d/ccres.doc",
    "https://usuaris.tinet.cat/debb/HGT/bmelic1.d/bmelic1.doc",
    "https://usuaris.tinet.cat/debb/HGT/bmelic2.d/bmelic2.doc",
    "https://usuaris.tinet.cat/debb/HGT/mloti.d/mloti.doc",
    "https://usuaris.tinet.cat/debb/HGT/smel.d/smel.doc",
    "https://usuaris.tinet.cat/debb/HGT/atum1c1.d/atum1c1.doc",
    "https://usuaris.tinet.cat/debb/HGT/atum1c2.d/atum1c2.doc",
    "https://usuaris.tinet.cat/debb/HGT/atum2c1.d/atum2c1.doc",
    "https://usuaris.tinet.cat/debb/HGT/atum2c2.d/atum2c2.doc",
    "https://usuaris.tinet.cat/debb/HGT/rpxx.d/rpxx.doc",
    "https://usuaris.tinet.cat/debb/HGT/rconorii.d/rconorii.doc",
    "https://usuaris.tinet.cat/debb/HGT/nmen1.d/nmen1.doc",
    "https://usuaris.tinet.cat/debb/HGT/nmen2.d/nmen2.doc",
    "https://usuaris.tinet.cat/debb/HGT/rsola.d/rsola.doc",
    "https://usuaris.tinet.cat/debb/HGT/rsolac2.d/rsolac2.doc",
    "https://usuaris.tinet.cat/debb/HGT/cjen.d/cjen.doc",
    "https://usuaris.tinet.cat/debb/HGT/hpyl.d/hpyl.doc",
    "https://usuaris.tinet.cat/debb/HGT/hpyl99.d/hpyl99.doc",
    "https://usuaris.tinet.cat/debb/HGT/baphi.d/baphi.doc",
    "https://usuaris.tinet.cat/debb/HGT/bsp.d/bsp.doc",
    "https://usuaris.tinet.cat/debb/HGT/ecoli.d/ecoli.doc",
    "https://usuaris.tinet.cat/debb/HGT/ecoli2.d/ecoli2.doc",
    "https://usuaris.tinet.cat/debb/HGT/ecoli3.d/ecoli3.doc",
    "https://usuaris.tinet.cat/debb/HGT/sent.d/sent.doc",
    "https://usuaris.tinet.cat/debb/HGT/styp.d/styp.doc",
    "https://usuaris.tinet.cat/debb/HGT/ypestis.d/ypestis.doc",
    "https://usuaris.tinet.cat/debb/HGT/hinf.d/hinf.doc",
    "https://usuaris.tinet.cat/debb/HGT/pmul.d/pmul.doc",
    "https://usuaris.tinet.cat/debb/HGT/paer.d/paer.doc",
    "https://usuaris.tinet.cat/debb/HGT/vcolc1.d/vcolc1.doc",
    "https://usuaris.tinet.cat/debb/HGT/vcolc2.d/vcolc2.doc",
    "https://usuaris.tinet.cat/debb/HGT/vvul1c1.d/vvul1c1.doc",
    "https://usuaris.tinet.cat/debb/HGT/vvul1c2.d/vvul1c2.doc",
    "https://usuaris.tinet.cat/debb/HGT/vvul2c1.d/vvul2c1.doc",
    "https://usuaris.tinet.cat/debb/HGT/vvul2c2.d/vvul2c2.doc",
    "https://usuaris.tinet.cat/debb/HGT/xcamp.d/xcamp.doc",
    "https://usuaris.tinet.cat/debb/HGT/xcitri.d/xcitri.doc",
    "https://usuaris.tinet.cat/debb/HGT/xfas.d/xfas.doc",
    "https://usuaris.tinet.cat/debb/HGT/bbur.d/bbur.doc",
    "https://usuaris.tinet.cat/debb/HGT/tpal.d/tpal.doc",
    "https://usuaris.tinet.cat/debb/HGT/tmar.d/tmar.doc",
    "https://usuaris.tinet.cat/debb/HGT/dra1.d/dra1.doc",
    "https://usuaris.tinet.cat/debb/HGT/dra2.d/dra2.doc"
]




def download_files(list_of_files):
    for link in list_of_files:
        #print(link)
        file_name = os.path.basename(link)
        file_name = file_name.replace(".doc", ".tsv")
        output_file = os.path.join(ROOT_FOLDER,file_name)
        #print(output_file)
        if os.path.isfile(output_file):
            print(f"{output_file} is downloaded")
        else:
            subprocess.run(['wget','-c', link ,'-O', output_file,], check=True)









if __name__ == '__main__':
    if not os.path.isdir(ROOT_FOLDER):
        os.makedirs(ROOT_FOLDER)
        print("created folder : ", ROOT_FOLDER)
    else:
        print(ROOT_FOLDER, "folder already exists.")
        
    download_files(LINK_TO_FILES)
        
    
    