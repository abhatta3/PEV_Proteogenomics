#https://pythonhosted.org/pyteomics/data.html
def make_list_mzid(inputevent,scanf,d_pep=1,d_index=2,d_fname=0,i_ev=0,i_seq=2,i_type=1):
    output=[]
    pep={}
    event_type={}
    novel={}
    for line in inputevent:
        if '#' in line:
            continue
        line=line.strip().split('\t')
        if len(line)<3:
            continue
        if '.' in line[i_seq]:
            pepseq=line[i_seq].split('.')[1].translate(None, '1234567890._:-*!@#$?')
        else:
            pepseq=line[i_seq].translate(None, '1234567890._:-*!@#$?')
            
        novel[pepseq]=line[i_seq]
        if pep.has_key(pepseq):
            old=pep[pepseq]
            old.append(line[i_ev])
            pep[pepseq]=old
        else:
            pep[pepseq]=[line[i_ev]]
        event_type[pepseq]=line[i_type]
    output.append("#Event\tIndex_number\tFile\tPeptide\tEvent_Type\n") 
           
    with open(scanf,"r") as sf:
        for line in sf:
            if 'PSM' not in line:
                continue
            line=line.strip().split('\t')
            if len(line)<8:
                continue
            if '.' in line[d_pep]:
                pepseq=line[d_pep].split('.')[1].translate(None, '1234567890._:-*!@#$?')
            else:
                pepseq=line[d_pep].translate(None, '1234567890._:-*!@#$?')
            fname=line[d_fname]
            if '=' in line[d_index]:
                index=line[d_index].split('=')[1]
            else:
                index=line[d_index]
                
            if pep.has_key(pepseq):
                pev=pep[pepseq]
                pis=pev[0]
                if len(pev)>1:
                    for itm in pev[1:]:
                        pis+=','+itm
                output.append("%s\t%s\t%s\t%s\t%s\n" %(pis,index,fname,line[d_pep],event_type[pepseq]))
    return output

