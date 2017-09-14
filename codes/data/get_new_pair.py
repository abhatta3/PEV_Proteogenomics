dic_ms={}
outf=open("pair.txt","w")
with open("inf1","r") as in1:
    for line in in1:
        if "#" in line:
            continue
        line=line.strip().split("\t")
        dic_ms[line[0]]=int(line[1])
    for w in dic_ms:
        wmass=dic_ms[w]
        for ww in dic_ms:
            if w != ww:
                outf.write("%s\t%s\t%d\n" %(w,ww,wmass-dic_ms[ww]))
        
