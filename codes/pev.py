import os, sys, getopt
import time
from utility import getannotdictionary
from novel_filtering import Filter
from Make_list import make_list
from extract_spec import extracter
from get_keywords import get_key
from modification_mutation_search import generate_alternatives
from read_mzid import mzid_tsv_convert
from Rescore import rescore

def print_help():
    print 'python pev.py -i <input event file> -s <input spectra files>  -d <input database search files> -o <outputfile>'

def main(argv):
    input_event = ''
    input_spectra = ''
    input_psm = ''
    try:
        opts, args = getopt.getopt(argv,"hi:s:d:o:")
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help()
            sys.exit()
        elif opt in "-i":
            input_event = arg
        elif opt in "-o":
            outfile = arg
        elif opt in "-d":
            input_psm = arg
        elif opt in "-s":
            input_spectra = arg
    if input_event == '' or outfile=='' or input_psm=='' or input_spectra=='':
        print("Error: Wrong parameter. File name missing") 
        print_help()
        sys.exit()
    cpath=os.path.dirname(os.path.realpath(sys.argv[0]))
    reffilename=cpath+"/data/alternate.fasta.tab"
    rmssfilename=cpath+"/data/rsiduemass"

    par_map_mgf={}

    if os.path.exists('params/params.xml'):
        with open('params/params.xml',"r") as param: 
            for line in param:
                if '"upload_file_mapping">mgffile' in line:
                    line=line.strip().split(">")[1]
                    new_mgf=line.split("|")[0]
                    old_mgf=line.split("/")[-2].split("<")[0]
                    par_map_mgf[old_mgf]=new_mgf

    ms_dict=getannotdictionary(rmssfilename)

    (novel_list,ref_list,itraq_tag)=Filter(input_event,reffilename)

    print("Novel filtering done.")
    make_list_psm_list=make_list(par_map_mgf,novel_list, input_psm, d_pep=8,d_index=1,d_fname=0,i_ev=0,i_seq=2,i_type=1)
    
    print("PSM list done")    
    singlemgf=outfile+"/single.mgf"
    
    extracted_psm_list=extracter(make_list_psm_list,input_spectra,singlemgf)

    print("single mgf done\n")

    start_time = time.time()
    
    key_out=get_key(extracted_psm_list,singlemgf,ms_dict,itraq_tag,4)

    print("Keyword Executes in (seconds) ---\t%s\n" %(time.time() - start_time))
    
    #print("keyword done\n")    
    
    singlefasta=outfile+"/single.fa"
    
    alt_outf=outfile+"/Alternatives.txt"

    start_time = time.time()
    
    generate_alternatives(key_out,reffilename,singlefasta,alt_outf,cpath)

    print("Alternative Executes in (seconds) ---\t%s\n" %(time.time() - start_time))

    #print("alternative search done")

    mzidout=outfile+"/msgfout.mzid"

    start_time = time.time()
    
    if itraq_tag==0:
        os.system("java -Xmx3500M -jar "+cpath+"/msgfplus/MSGFPlus.jar  -d  "+singlefasta+" -s  "+singlemgf+"  -o "+mzidout+" -t 1000ppm -m 0 -inst 0 -e 1 -ti -1,2 -ntt 2 -tda 0 -n 8 -thread 7  -mod  "+cpath+"/data/Modifications_msgf.txt")
    elif itraq_tag==4:
        os.system("java -Xmx3500M -jar "+cpath+"/msgfplus/MSGFPlus.jar  -d  "+singlefasta+" -s  "+singlemgf+"  -o "+mzidout+" -t 1000ppm -m 0 -inst 0 -e 1 -ti -1,2 -ntt 2 -tda 0 -n 8 -thread 7  -mod  "+cpath+"/data/Modifications_msgf_itq4.txt")
    elif itraq_tag==8:
        os.system("java -Xmx3500M -jar "+cpath+"/msgfplus/MSGFPlus.jar  -d  "+singlefasta+" -s  "+singlemgf+"  -o "+mzidout+" -t 1000ppm -m 0 -inst 0 -e 1 -ti -1,2 -ntt 2 -tda 0 -n 8 -thread 7  -mod  "+cpath+"/data/Modifications_msgf_itq8.txt")

    print("MSGF+ search Executes in (seconds) ---\t%s\n" %(time.time() - start_time))

    
    #print("MSGF search done")
    tsvout=outfile+"/msgfout.tsv"
    mzid_tsv_convert(mzidout,tsvout)
    print("mzid to tsv done")
    outputfile=outfile+"/result.txt"

    start_time = time.time()
    
    rescore(tsvout,outputfile,extracted_psm_list,singlemgf,novel_list,key_out,itraq_tag,ref_list) #,pep_by_kw)

    print("Rescore Executes in (seconds) ---\t%s\n" %(time.time() - start_time))

    
   
            
if __name__ == "__main__":
    start_time = time.time()
    main(sys.argv[1:])
    print("Executes in --- %s seconds ---" %(time.time() - start_time))
