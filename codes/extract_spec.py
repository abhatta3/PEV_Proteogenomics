#https://pythonhosted.org/pyteomics/data.html
import os
import os.path

def extracter(inputf,spec_dir,spec_outfile):
    soutput=open(spec_outfile,"w")
    output=[]
    event_scan={}
    for line in inputf:
        if '#' in line:
            continue
        line=line.strip().split('\t')
        if len(line)>=5:
            index=int(line[1])
            fname=line[2]
            pev=line[0]
            peptide=line[3]
            etype=line[4]
            if event_scan.has_key(fname):
                sdic=event_scan[fname]
                sdic[index]={'pev':pev,'pseq':peptide,'etype':etype}
                event_scan[fname]=sdic
            else:
                temp={'pev':pev,'pseq':peptide,'etype':etype}
                sdic={index:temp}
                event_scan[fname]=sdic
    newscanno=0
    output.append("#Event\tIndex_number\tFile\tPeptide\tEvent_Type\n")            
    for fname in event_scan:
        scanlist=event_scan[fname]
        #PATH='mgffile/'+fname
        if os.path.isdir(spec_dir):
            PATH=spec_dir+'/'+fname
            print "AS DIRECTORY"
        else:
            #PATH=os.path.dirname(os.path.abspath(spec_dir))
            PATH=fname
            print "AS PATH"

        flag=0
        if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
            with open(PATH,"r") as mgf:
                scanindex=0
                for line in mgf:
                    if len(line)>0:
                        if 'BEGIN IONS' in line:
                            if scanlist.has_key(scanindex):
                                pev=event_scan[fname][scanindex]['pev']
                                output.append("%s\t%d\t%s\t%s\t%s\n" %(pev,newscanno,spec_outfile,event_scan[fname][scanindex]['pseq'],event_scan[fname][scanindex]['etype']))
                                newscanno+=1
                                soutput.write(line)
                                flag=1
                            scanindex+=1
                        elif flag==1 and 'END IONS' not in line:
                            soutput.write(line)
                        elif flag==1 and 'END IONS' in line:
                            soutput.write(line)
                            flag=0
    soutput.close()
    return output
                             
