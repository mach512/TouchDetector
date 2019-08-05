import h5py
import numpy as np
from neuron import h, gui
from SampPurkinje import Purkinje
purkcell=Purkinje()
f=h5py.File('scaffold_full_IO_200.0x200.0_v3.hdf5','r')
cellpositions=np.array(f['positions'])
pcpositions=cellpositions[np.where(cellpositions[:,1]==4.0)]
h.define_shape()
pcsecall=h.SectionList()
pcsecall.wholetree(purkcell.soma)
listgr=[]
sectiondct={}
maincelldct=[]
dictcnt=0

for ij in range(len(pcpositions)):
	#exec('sectdict_%d={}'%ij)
	sectiondct={}
	for sec in pcsecall: 
		segstr=str(sec)
		sec.push()
		listgrX=[]
		listgrY=[]
		listgrZ=[]
		for i in range(int(h.n3d())): 
			listgrX.append(h.x3d(i))
			listgrY.append(h.y3d(i))
			listgrZ.append(h.z3d(i))
		
		if str(sec)=="axonAIS" or str(sec)=="axonmyelin" or str(sec)=="axonmyelin4" or str(sec)=="axonNOR3" or str(sec)=="axoncoll2" or str(sec)=="axonmyelin3" or str(sec)=="axoncoll" or str(sec)=="axonNOR2" or str(sec)=="axonmyelin2" or str(sec)=="axonNOR" or str(sec)=="axonAISK" or str(sec)=="axonAIS":
			segstr=str(sec)
			buf=segstr
		else: 
			segstr=segstr.split('.')
			segstr=segstr[2].split('[')
			if segstr[0]=='soma': 
				buf=segstr[0]
			else: 
				secname=segstr[0]
				segstr=segstr[1].split(']')
				secid=segstr[0]
				buf=secname+"_"+secid
		#print 'buf contains:', buf
		
		maxX=max(listgrX)+pcpositions[ij,2]
		minX=min(listgrX)+pcpositions[ij,2]
		maxY=max(listgrY)+pcpositions[ij,3]
		minY=min(listgrY)+pcpositions[ij,3]
		maxZ=max(listgrZ)+pcpositions[ij,4]
		minZ=min(listgrZ)+pcpositions[ij,4]
		
		#getallx=np.column_stack((minX, maxX))
		#getally=np.column_stack((minY, maxY))
		#getallz=np.column_stack((minZ, maxZ))
		
		#h1=h5py.File('newpurkdata.hdf5', 'a')
		#h1.create_dataset('xcoord', data=getallx)
		#h1.create_dataset('ycoord', data=getally)
		#h1.create_dataset('zcoord', data=getallz)
		
		
		with open('savepcsegmentcoords.dat', 'a') as fsavaxon:
			fsavaxon.write(str(ij))
			fsavaxon.write('\t')
			fsavaxon.write(str(buf))
			fsavaxon.write('\t')
			fsavaxon.write(str(minX))
			fsavaxon.write('\t')
			fsavaxon.write(str(minY))
			fsavaxon.write('\t')
			fsavaxon.write(str(minZ))
			fsavaxon.write('\n')
			fsavaxon.write(str(ij))
			fsavaxon.write('\t')
			fsavaxon.write(str(buf))
			fsavaxon.write('\t')
			fsavaxon.write(str(maxX))
			fsavaxon.write('\t')
			fsavaxon.write(str(maxY))
			fsavaxon.write('\t')
			fsavaxon.write(str(maxZ))
			fsavaxon.write('\n')
		h.pop_section()
#h1.close()