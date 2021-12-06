#!/usr/bin/env python3
d={}
file=open("Probes.txt","r")
rows=file.readlines()
for i in range(0,len(rows)):
	row=rows[i].rstrip()
	row=row.split("\t")
	id=row[0]
	if id in d:
		continue
	else:
		d[id]=row[1]
		#d1[id]=[[row1[1],row1[2]]]

#print(d)
#UpregulatedinHCTvsHEK293T.txt

file=open("UpregulatedinH9vsHEK293T.txt","r")
rows=file.readlines()
for i in range(0,len(rows)):
        row=rows[i].rstrip()
        row=row.split("\t")
        id2=row[0]
        if id2 in d:
                print(row[0]+"\t"+d[id2]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3]+"\t"+row[4]+"\t"+row[5]+"\t"+row[6])
  
