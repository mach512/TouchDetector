from __future__ import division
from ipyparallel import Client
import pandas as pd
import json
from mpi4py import MPI
import h5py
import numpy as np
rc=Client(profile='setupparallelcluster')
dview=rc[:]
dview.execute('from __future__ import division')
dview.execute('import numpy as np')
dview.execute('import pandas as pd')
dview.execute('import h5py')
dview.execute('import matplotlib.pyplot as plt')
dview.execute('import pandas as pd')
dfgrc=pd.read_csv('singlegrcsecdata.dat',sep='\t',names=["gid","xdata","ydata","zdata"])
dfgoc=pd.read_csv('singlegocsecdata.dat',sep='\t',names=["gid","xdata","ydata","zdata"])
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()
#with dview.sync_imports(): 
#	import pandas as pd
#	dfgoc=pd.read_csv('singlegocsecdata.dat',sep='\t',names=["gid","xdata","ydata","zdata"])

f=h5py.File('scaffold_full_dcn_400.0x400.0_v3.hdf5','r')
cellpositions=np.array(f['positions'])
grcpositions=cellpositions[np.where(cellpositions[:,1]==3.0)]
gocpositions=cellpositions[np.where(cellpositions[:,1]==1.0)]
dictgrcXstore={}
dictgrcYstore={}
dictgrcZstore={}
dictgrcgidstore={}
for ii in range(5): 
	newdfgrcX=dfgrc["xdata"]+grcpositions[ii,2]
	newdfgrcY=dfgrc["ydata"]+grcpositions[ii,4]
	newdfgrcZ=dfgrc["zdata"]+grcpositions[ii,3]
	dictgrcgidstore[ii]=dfgrc["gid"]
	dictgrcXstore[ii]=newdfgrcX
	dictgrcYstore[ii]=newdfgrcY
	dictgrcZstore[ii]=newdfgrcZ
#dictallstoregrc=[dictgrcXstore[0],dictgrcYstore[0],dictgrcZstore[0]]
dictgocXstore=dfgoc["xdata"]+gocpositions[0,2]
dictgocYstore=dfgoc["ydata"]+gocpositions[0,4]
dictgocZstore=dfgoc["zdata"]+gocpositions[0,3]
newdictgocxstor={}
newdictgocystor={}
newdictgoczstor={}
dictgocgidstore={}
newdictgocxstor[0]=dictgocXstore
newdictgocystor[0]=dictgocYstore
newdictgoczstor[0]=dictgocZstore
dictgocgidstore[0]=dfgoc["gid"]
print 'length of dictgocgidstore[0] up above:', len(dictgocgidstore[0])
#print newdictgocxstor
dview.push(newdictgocxstor)
dview.push(newdictgocystor)
dview.push(newdictgoczstor)
dview.push(dictgocgidstore)
print 'next step'
checkdict={'43':np.array([1,2,3]),'34':np.array([4,5,6]),'343':np.array([7,8,9])}
#this is for the scattering dictionaries
def scatter_dict(view,name,d): 
	ntargets=len(view)
	keys=d.keys()
	for i, target in enumerate(view.targets): 
		subd={}
		for key in keys[i::ntargets]:
			subd[key]=d[key]
		view.client[target][name]=subd

scatter_dict(dview,'test',dictgrcXstore)
scatter_dict(dview,'test1',dictgrcYstore)
scatter_dict(dview,'test2',dictgrcZstore)
scatter_dict(dview,'gidtest',dictgrcgidstore)
#print 'size of dictgrcXstore:', len(dictgrcXstore),len(dictgrcXstore[0])
def gather_gidgrc_dict(view,name): 
	merged={}
	cellgidtotal=[]
	for d in view.pull(name): 
		listgidgrc=list(d.values())
		listgidgoc=list(dictgocgidstore.values())
	#	print 'length of the listgidgrc:', len(listgidgrc)
		if len(listgidgrc)>0:
			cellgiddata=[(x1,y1) for x,x1 in enumerate(listgidgrc[0]) for y,y1 in enumerate(listgidgoc[0])]
			cellgidtotal.append(cellgiddata)
		merged.update(d)
	return merged,cellgidtotal
resgid,getgiddata=gather_gidgrc_dict(dview,'gidtest')
def gather_dict(view,name): 
	merged={}
	cellxtotal=[]
	cellnumxtotal=[]
	cnt=0
	for d in view.pull(name):
		#cellxdata=[(x1,y1) for x,x1 in d.items() for y,y1 in newdictgocxstor.items()]
		#for key in d.keys(): 
		listdgrc=list(d.values())
		listdgoc=list(newdictgocxstor.values())
		if len(listdgrc) > 0:
			cellxdata=[(x,y,x1,y1) for x,x1 in enumerate(listdgrc[0]) for y,y1 in enumerate(listdgoc[0])]
			cellxtotal.append(cellxdata)
		cnt+=1
		merged.update(d)
	return merged,cellxtotal
def gather_dict_Y(view,name): 
	merged={}
	cellytotal=[]
	cnt=0
	for d in view.pull(name):
		listdgrc=list(d.values())
		listdgoc=list(newdictgocystor.values())
		if len(listdgrc) > 0:  
			cellydata=[(x1,y1) for x,x1 in enumerate(listdgrc[0]) for y,y1 in enumerate(listdgoc[0])]
			cellytotal.append(cellydata)
		cnt+=1
		merged.update(d)
	return merged,cellytotal
def gather_dict_Z(view,name): 
	merged={}
	cellztotal=[]
	cnt=0
	for d in view.pull(name):
		listdgrc=list(d.values())
		listdgoc=list(newdictgoczstor.values())
		if len(listdgrc) > 0: 
			cellzdata=[(x1,y1) for x,x1 in enumerate(listdgrc[0]) for y,y1 in enumerate(listdgoc[0])]
			cellztotal.append(cellzdata)
		cnt+=1
		merged.update(d)
	return merged,cellztotal
res,getxdata=gather_dict(dview,'test')
resY,getydata=gather_dict_Y(dview,'test1')
resZ,getzdata=gather_dict_Z(dview,'test2')
#@dview.parallel(block=True)
def calcintersection(a,b,c,d,cellid,grcgid,gocgid): 
	marginofError=0.1
	r=np.subtract(b,a)
	s=np.subtract(d,c)
	q=np.subtract(a,c)
	dotqr=np.dot(q,r)
	dotqs=np.dot(q,s)
	dotrs=np.dot(r,s)
	dotrr=np.dot(r,r)
	dotss=np.dot(s,s)
	denom=dotrr*dotss-dotrs*dotrs
	numer=dotqs*dotrs-dotqr*dotss
	#print 'denom:', denom
	if denom!=0: 
		t=numer/denom
		u=(dotqs+t*dotrs)/dotss
		p0=a+t*r
		p1=c+u*s
		onsegment=False
		intersects=False
		mask=np.logical_and(t>=0,t<=1)
		mask1=np.logical_and(u>=0,u<=1)
		if mask==True and mask1==True: 
			#print 'mask and mask1 are true'
			onSegment=True
		res=p0-p1
		mag=np.sqrt(res.dot(res))
		if (mag < marginofError): 
			#print 'intersects is true'
			intersects=True
			with open('checkingtrueconnections.dat','a') as fsave: 
				fsave.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(p0,p1,mag,cellid,grcgid,gocgid))
		if onsegment==True and intersects==True: 
			print 'both are true'
##############################
cnt=0
cnt1=-1
mocnt=0

for ij in range(len(getxdata)):
	print 'length of getxdata:', len(getxdata),len(getxdata[0]),getxdata[0][0][3]
	mocnt=0
	for jk in range(len(getxdata[ij])-1): 
		if cnt==0:
			storeval=getxdata[ij][jk][3]
			cnt+=1
		if jk > 0:
			if getxdata[ij][jk][3]==storeval: 
				print 'jk value here:', jk
				grcstart=[getxdata[ij][jk-1][2],getydata[ij][jk-1][0],getzdata[ij][jk-1][0]]
				grcstop=[getxdata[ij][jk][2],getydata[ij][jk][0],getzdata[ij][jk][0]]
				print 'grcstart',grcstart,'grcstop',grcstop
				itemgrcfind=getgiddata[ij][0][mocnt]
				getgidnum=[getgiddata[ij][jk][0]]
				#print 'getgidnum:', getgidnum
				#if grcstart==grcstop:
					#print grcstart,grcstop
				mocnt+=1
				cnt1=0
			if cnt1==0: 
				if grcstart!=grcstop:
					#if grcstart[0] in dictgrcXstore:
					#next(item for item in dictgrcXstore.values() if item==grcstart[0])
			#newlistgrc=list(dictgrcXstore.values())
					#if grcstart[0] in newlistgrc: 
						#getitemgrc=[item for item in newlistgrc[0] if item==grcstart[0]]
					#	print grcstart[0],len(newlistgrc[4])
					#for ij in range(len(newlistgrc)):
					#	getitemgrc=[(i,x) for i, x in enumerate(newlistgrc[0]) if x==grcstart[0]]
					#	print getitemgrc
					#if grcstart[0] in newlistgrc[0]: 
					#	print 'grcstart[0] is:', grcstart[0]
					#	print 'length of dictgrcXstore:', newlistgrc[0].index[grcstart[0]]
					gocstart=[getxdata[ij][jk][3],getydata[ij][jk][1],getzdata[ij][jk][1]]
					gocstop=[getxdata[ij][jk+1][3],getydata[ij][jk+1][1],getzdata[ij][jk+1][1]]
					#print 'gocstart:', gocstart,'gocstop:', gocstop
					calcintersection(grcstart,grcstop,gocstart,gocstop,ij,getxdata[ij][jk][0],getxdata[ij][jk][1])
				#print 'jk value:', jk
				#print 'grcstart:',grcstart,'grcstop:',grcstop,'gocstart:',gocstart,'gocstop:',gocstop
					#cnt1+=1
		# if ij%2==0:
			# grcstart=[getxdata[ij][jk][0],getydata[ij][jk][0],getzdata[ij][jk][0]]
			# grcstart=[res[0][ij],resY[0][ij],resZ[0][ij]]
			# grcstop=[res[0][ij+1],resY[0][ij+1],resZ[0][ij+1]]
			# gocstart=[getxdata[ij][jk][1],getydata[ij][jk][1],getzdata[ij][jk][1]]
			# gocstop=[getxdata[ij][jk+1][1],getydata[ij][jk+1][1],getzdata[ij][jk+1][1]]
			# grcstop=[getxdata[ij][jk+1][0],getydata[ij][jk+1][0],getzdata[ij][jk+1][0]]
			#calcintersection(grcstart,grcstop,gocstart,gocstop)
		
#print 'grcstart:', getxdata[0][176][0],getydata[0][176][0],getzdata[0][176][0]
#print 'grcstop:', getxdata[0][353][0],getydata[0][353][0],getzdata[0][353][0]
#print 'gocstart:', getxdata[0][0][1],getydata[0][0][1],getzdata[0][0][1]
#print 'gocstop:', getxdata[0][1][1],getydata[0][1][1],getzdata[0][1][1]
#print 'grcstop:', rcstop
#dview.scatter('dictgrcXStore',dictgrcXstore)
#dview.scatter('dictgrcYStore',dictgrcYstore)
#dview.scatter('dictgrcZStore',dictgrcZstore)
print 'distributed the dictgrcXstore'
