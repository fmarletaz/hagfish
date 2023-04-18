 #!/usr/bin/env python

import sys,csv
import pyfastx
from collections import defaultdict
from collections import Counter

clgAlis=defaultdict(list)


for r in csv.reader(open(sys.argv[1]),delimiter='\t'):
    if r[0]=='fid': continue
    fid,clg,bfl,spu,bla,sko,pfl,apl=r[0:8]
    #clgAlis[clg]=defaultdict(list)
    select={}
    for spn,rsp in zip(['Braflo','Strpur','Bralan','Sackow','Ptyfla','Acaplan'],[bfl,spu,bla,sko,pfl,apl]):
        if not rsp=='':
            select[spn+'_'+rsp]='NA'
    #select={'Braflo_'+bfl:'NA'}
    #print(select)
    for c in r[8:]:
        #print(len(c))
        if len(c)==0: continue
        for g in c.split(';'):
            gid,cld=g.split('|')
            #print(gid,cld)
            select[gid]=cld
    ffid=fid.split('.')[0]
    alFaDir="/Users/fmarletaz/Dropbox/Genomes/Hagfish/Gene_families/aliFasta2"
    #alFaDir="/Users/fmarletaz/Desktop/aliFasta2"
    for name, seq in pyfastx.Fasta("{}/{}.al.cl.fa".format(alFaDir,ffid), build_index=False):    
        if name in select:
            #print(name,select[name])
            sp=name.split('_',1)[0]
            txn="{}_{}".format(sp,select[name])
            clgAlis[clg].append((fid,clg,txn,seq))
            #clgAlis[clg][txn].append(seq)

for clg,alis in clgAlis.items():
    Txns=set([c[2] for c in alis])
    Alxs=defaultdict(dict)
    for a in alis:
        Alxs[a[0]][a[2]]=a[3]
    print(clg,len(Txns),len(Alxs))
    Sqxs=defaultdict(str)
    Ctxs=defaultdict(int)
    for fid,alx in Alxs.items():
        if len(alx)/len(Txns)<=0.2:continue
        #print(clg,fid,len(alx),len(Txns))

        for tx in Txns:
            if tx in Alxs[fid]:
                #print(tx)
                Sqxs[tx]+=Alxs[fid][tx]
                Ctxs[tx]+=1
            else:
                #print('no',tx)
                Sqxs[tx]+='X' * len(list(Alxs[fid].values())[0])
                #print('X' * len(list(Alxs[fid].values())[0]))

    with open('{}.fa'.format(clg),'w') as out:
        for tx in Sqxs:
            #print(tx,Ctxs[tx],len(Txns))
            if Ctxs[tx]<0.001*len(Txns) or len(list(Counter(Sqxs[tx]).keys()))<5:
                print('   excluding',tx)
                continue
            #print('  ',tx,Ctxs[tx],)
            if not len(list(Counter(Sqxs[tx]).keys()))<5:
                out.write(">{}\n{}\n".format(tx,Sqxs[tx]))
            else:
                print('  >',fid,tx,'skipped!')


#for clg in clgAlis: 
#    with open('{}.fa'.format(clg),'w') as out:
#        for txn in clgAlis[clg]:
