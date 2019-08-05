#implementation of the algorithm. 
#loading the scaffold files (hdf5) to get the positions for moving the Purkinje neuron to the appropriate position. 
import pandas as pd
import h5py
import os
import timeit
from neuron import h, gui
from SampPurkinje import Purkinje
import multiprocessing as mp
import numpy as np
fid='GrcTouchLocations.dat'
from neuron.rxd.morphology import parent, parent_loc
import json
pulist=[]
h.load_file('New_SMGrC.hoc')
h.load_file('SM_golgi.hoc')

#st=pd.HDFStore('data.h5',mode='w')

#from record3dseglocs import GetNeuronData
grcgid=[]
pcgid=[]

f=h5py.File('scaffold_full_IO_200.0x200.0_v3.hdf5', 'r')
cellpositions=np.array(f['positions'])
gocpositions=cellpositions[np.where(cellpositions[:,1]==1.0)]
grcpositions=cellpositions[np.where(cellpositions[:,1]==3.0)]
pcpositions=cellpositions[np.where(cellpositions[:,1]==4.0)]

length_range=len(grcpositions)
getxsmin=[grcpositions[i,2]+(-0.21148) for i in range(length_range)]
getxsmax=[grcpositions[i,2]+0.139 for i in range(length_range)]
getysmin=[grcpositions[i,4]+(-2.8867) for i in range(length_range)]
getysmax=[grcpositions[i,4]+2.601 for i in range(length_range)]
getzsmin=[grcpositions[i,3]+0 for i in range(length_range)]
getzsmax=[grcpositions[i,3]+0 for i in range(length_range)]

#print 'length of the getxsmin:', len(getxsmin)
#print 'length of grcpositions:', len(grcpositions)
#print 'length of pcposition:', len(pcpositions)

#Function to identify sections and change their positions. 
def setup_grcmorphology(i,pursections,x,y,z):
	if i > 3500: 
		print 'i value is:', i
	cell=h.SMGrC2019()
	print 'grc here'
	h.define_shape()
	grcallsec=h.SectionList()
	grcallsec.wholetree(cell.soma[0])
	cell.soma[0].push()
	n=h.n3d(sec=cell.soma[0])
	xs=[h.x3d(i) for i in range(int(n))]
	ys=[h.y3d(i) for i in range(int(n))]
	zs=[h.z3d(i) for i in range(int(n))]
	ds=[h.diam3d(i) for i in range(int(n))]
	j=0
	#sec.push()
	for a, b, c, d in zip(xs,ys,zs,ds):
		#print 'sec here is:', sec 
		h.pt3dchange(j,a+x,b+y,c+z,d)
		j+=1 
	h.define_shape()
	h.pop_section()
	for sec in pursections:
		xsecappend=[]
		ysecappend=[]
		zsecappend=[]
		segstr=str(sec)
		#print 'segstr:', segstr
		segstr=segstr.split('.')
		if str(sec)=="axonAIS" or str(sec)=="axonmyelin" or str(sec)=="axonmyelin4" or str(sec)=="axonNOR3" or str(sec)=="axoncoll2" or str(sec)=="axonmyelin3" or str(sec)=="axoncoll" or str(sec)=="axonNOR2" or str(sec)=="axonmyelin2" or str(sec)=="axonNOR" or str(sec)=="axonAISK" or str(sec)=="axonAIS":
			segstr=str(sec)
			ext="axon"
			segstr=segstr[:segstr.find(ext)+len(ext)]
			#print 'the section name here is axonAIS'
		#print 'before segstr[2] split', segstr
		else: 
			segstr=segstr[2].split('[')
			segstr=''.join([i for i in segstr[0] if not i.isdigit()])
		if segstr=='bs': 
			#print 'this section is dendrite', sec
			sec.push()
			for kj in range(int(h.n3d()-1)): 
				xsecappend.append(h.x3d(kj))
				#append_secX_list.append(h.x3d(h.n3d()-1))
				ysecappend.append(h.y3d(kj))
				#append_secY_list.append(h.y3d(h.n3d()-1))
				zsecappend.append(h.z3d(kj))
				purkallcoord=[sec.x3d(kj),sec.y3d(kj),sec.z3d(kj)]
				purkradius=(sec.diam/2)
				for gsec in grcallsec: 
					gsegstr=str(gsec)
					gsegstr=gsegstr.split('.')
					gsegstr=gsegstr[1].split('[')
					gsegid=gsegstr[1].split(']')
					if gsegstr[0]=='axon': 
						if gsegid > 1 and gsegid < 20: 
							for jj in range(int(gsec.n3d()-1)): 
								gxmin=gsec.x3d(jj)
								gymin=gsec.y3d(jj)
								gzmin=gsec.z3d(jj)
								grallcoord=[gsec.x3d(jj),gsec.y3d(jj),gsec.z3d(jj)]
								dst=distance.euclidean(grallcoord,purkallcoord)
								granradius=(gsec.diam/2)
								if dst < (purkradius+granradius): 
									print 'probable connection'
								
				#append_secZ_list.append(h.z3d(h.n3d()-1))
			h.pop_section()
		
def callback(result): 
	print 'rsult:', result
	
def func_local_grcneurons(gid,sections): 
	print 'gid value is:', gid
	#getpurkneuron=pulist[gid]
	#print 'getpurkneuron:', getpurkneuron
	append_secX_list=[]
	append_secY_list=[]
	append_secZ_list=[]
	hPurk=h5py.File('dataPurk.hdf5','a')
	for sec in sections:
		segstr=str(sec)
		#print 'segstr:', segstr
		segstr=segstr.split('.')
		if str(sec)=="axonAIS" or str(sec)=="axonmyelin" or str(sec)=="axonmyelin4" or str(sec)=="axonNOR3" or str(sec)=="axoncoll2" or str(sec)=="axonmyelin3" or str(sec)=="axoncoll" or str(sec)=="axonNOR2" or str(sec)=="axonmyelin2" or str(sec)=="axonNOR" or str(sec)=="axonAISK" or str(sec)=="axonAIS":
			segstr=str(sec)
			ext="axon"
			segstr=segstr[:segstr.find(ext)+len(ext)]
			#print 'the section name here is axonAIS'
		#print 'before segstr[2] split', segstr
		else: 
			segstr=segstr[2].split('[')
			segstr=''.join([i for i in segstr[0] if not i.isdigit()])
		if segstr=='bs': 
			#print 'this section is dendrite', sec
			sec.push()
			for kj in range(int(h.n3d()-1)): 
				append_secX_list.append(h.x3d(kj))
				#append_secX_list.append(h.x3d(h.n3d()-1))
				append_secY_list.append(h.y3d(kj))
				#append_secY_list.append(h.y3d(h.n3d()-1))
				append_secZ_list.append(h.z3d(kj))
				
				#append_secZ_list.append(h.z3d(h.n3d()-1))
			h.pop_section()
	#find max and min from each list declared above
	purksecdata=np.column_stack((append_secX_list,append_secY_list,append_secZ_list))
	hPurkdset=hPurk.create_dataset('data_purk',data=purksecdata)
	print 'sec is:', sec
	#print 'first element:', append_secZ_list[0]
	min_secX_list=min(append_secX_list)
	max_secX_list=max(append_secX_list)
	min_secY_list=min(append_secY_list)
	max_secY_list=max(append_secY_list)
	min_secZ_list=min(append_secZ_list)
	max_secZ_list=max(append_secZ_list)
	
	#checking to see if the grcpositions x and y min, max fall within the boundaries of Purkinje neuron
	print 'min_secX_list,max_secX_list,min_secY_list,max_secY_list,min_secZ_list,max_secZ_list', min_secX_list,max_secX_list,min_secY_list,max_secY_list,min_secZ_list,max_secZ_list
	chkXcnt=0
	localitercnt=0
	grcxminlist=[]
	grcxmaxlist=[]
	grcyminlist=[]
	grcymaxlist=[]
	grczminlist=[]
	grczmaxlist=[]
	for (a,b,c,d,e,f) in zip(getxsmin, getxsmax,getysmin,getysmax,getzsmin,getzsmax):
		#print a,b,c,d
		localitercnt+=1
		if a > min_secX_list and b < max_secX_list and c > min_secY_list and d < max_secY_list and e > min_secZ_list and f < max_secZ_list:
			chkXcnt+=1 #no. of granule neurons that fall into boundary for each Purkinje neuron
			grcgid.append(localitercnt)
			pcgid.append(gid)
			grcxminlist.append(a)
			grcxmaxlist.append(b)
			grcyminlist.append(c)
			grcymaxlist.append(d)
			grczminlist.append(e)
			grczmaxlist.append(f)
		#if c > min_secY_list and d < max_secY_list: 
		#	chkXcnt+=1
	print 'length of grcgid:', len(grcgid)
	#for jj,iterjj in enumerate(grcgid): 
	#	print 'grcgid:', jj, iterjj
	print '21849 value:', grcpositions[0,2]
	data_grc={'grcgid':grcgid}
	data_pc={'pcgid':pcgid}
	data_xyz=np.column_stack((grcxminlist,grcxmaxlist,grcyminlist,grcymaxlist,grczminlist,grczmaxlist))
	#print 'data:', data
	h1=h5py.File('data.hdf5','a')
	data=np.column_stack((grcgid,pcgid))
	#print 'data:', data
	#df1=pd.DataFrame.from_dict(data_grc)#({'grcgid':grcgid,'pcgid':pcgid})
	#df2=pd.DataFrame.from_dict(data_pc)
	
	dset=h1.create_dataset('data_grc',data=data)
	dset1=h1.create_dataset('grcXYZ',data=data_xyz)
	h1.close()
	del grcxminlist
	del grcxmaxlist
	del grcyminlist
	del grcymaxlist
	del grczminlist
	del grczmaxlist
	del min_secX_list
	del max_secX_list
	del min_secY_list
	del max_secY_list
	del min_secZ_list
	del max_secZ_list
	del append_secX_list
	del append_secY_list
	del append_secZ_list
	#grcs=[pool.apply_async(setup_grcmorphology,args=(jj, sections, grcpositions[grcgid[jj],2], grcpositions[grcgid[jj],4], grcpositions[grcgid[jj],3])) for jj in range(len(grcgid)-1)]
	grcs=[pool.apply(callback,args=(i)) for i in range(10)]
	pool.close()
	print 'grcs:', grcs[:20]
	#grcs=[setup_grcmorphology(jj,sections,grcpositions[grcgid[jj],2], grcpositions[grcgid[jj],4], grcpositions[grcgid[jj],3]) for jj in range(len(grcgid)-1)]
	#print 'df1:', type(df1)
	#st=pd.HDFStore('data.h5',mode='a')
	#st.append('df',df1,data_columns=['grcgid'],index=False)
	#st.append('df',df2,data_columns=['pcgid'],index=False)
	#df1.to_hdf('data.h5', key='df', table=True, data_columns=['grcgid','pcgid'],mode='a')
	#st.close()
	print 'chkXcnt:', chkXcnt
	#print 'min granule x is:', getxsmin

def setup_pcmorphology(kk,x,y,z): 
	print 'i value is:', kk
	cellpc=Purkinje()
	h.define_shape()
	
	sections=h.SectionList()
	sections.wholetree(cellpc.soma)
	cellpc.soma.push()
	n=h.n3d(sec=cellpc.soma)
	xs=[h.x3d(i) for i in range(int(n))]
	ys=[h.y3d(i) for i in range(int(n))]
	zs=[h.z3d(i) for i in range(int(n))]
	ds=[h.diam3d(i) for i in range(int(n))]
	j=0
	#sec.push()
	for a, b, c, d in zip(xs,ys,zs,ds):
		#print 'sec here is:', sec 
		h.pt3dchange(j,a+x,b+y,c+z,d)
		j+=1 
	h.define_shape()
	h.pop_section()
	pulist.append(cellpc)
	#call another function to locate the local granule neurons
	cellpc.soma.push()
	getpurksecs=h.SectionList()
	getpurksecs.wholetree(cellpc.soma)
	func_local_grcneurons(kk,getpurksecs)

if __name__== '__main__': #this is extremely important for python scripts in windows
	filePath='C:/Users/mach5/OneDrive/Desktop/SM_granule_golgi_models/data.hdf5'
	file2='C:/Users/mach5/OneDrive/Desktop/SM_granule_golgi_models/dataPurk.hdf5'
	if os.path.exists(filePath): 
		os.remove(filePath)
	else: 
		print 'cannot delete file as it is not there'
	if os.path.exists(file2): 
		os.remove(file2)
	else: 
		print 'file 2 is not there'
		
	#del data.hdf5
	#del dataPurk.hdf5
	start=timeit.default_timer()
	pool=mp.Pool(8)
	#grcs=[pool.apply_async(setup_grcmorphology,args=(i, grcpositions[i,2], grcpositions[i,4], grcpositions[i,3])) for i in range(100)]
	#pcs=[pool.apply_async(setup_pcmorphology,args=(i,pcpositions[i,2],pcpositions[i,4],pcpositions[i,3])) for i in range(10)]
	
	pcs=[setup_pcmorphology(i,pcpositions[i,2],pcpositions[i,4],pcpositions[i,3]) for i in range(1)]
	stop=timeit.default_timer()
	print('Time:', stop-start)
	#pool.close()
	#pool.join()
	#grcs=[setup_grcmorphology(i,grcpositions[i,2],grcpositions[i,4],grcpositions[i,3]) for i in range(100)]
	slt=h.Shape()
	slt.view()
