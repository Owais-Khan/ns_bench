import matplotlib.pyplot as plt
import sys
import numpy as np
import operator
import sys
import glob

class CYLINDER_PLOT():
	def __init__(self):
		#Load all the cases
		self.cases=sys.argv[1:]
		
		#Get the tecplot files in each folder
		file_names=[glob.glob(i+"*.dat")[0] for i in sys.argv[1:]]
		
		#Open the files
		self.files=[open(i,'r') for i in file_names]
		self.label   =["ABI_CN_HR","ABI_CN_NR","ABI_CN_LR"]
		
		#Read headers in each file (first 3 lines)
		[i.readline() for i in self.files]
		[i.readline() for i in self.files]
		[i.readline() for i in self.files]

	def plot_data(self):
		#Create figures and axes for each case
		FIG=[]; AX=[]
		for i in range(0,1):
			Fig,ax=plt.subplots(1, sharex=True, sharey=False)
			FIG.append(Fig)
			AX.append(ax)

		#Loop over cases
		count_FILE=0
		count_file=0
		for FILE in self.files:
			#Loop over lines for each case
			t=[]; dp=[]; cd=[]; cl=[]
			for LINE in FILE:
				line=[float(i) for i in LINE.split()]
				t.append(line[0])
				dp.append(line[1])
				cd.append(-1*line[3])
				cl.append(-1*line[4])

			
			#Maximum and Minimum of Pressure, drag and lift.			
			dp_min=min(dp[int(0.5*len(t)):])
			dp_max_i, dp_max = max(enumerate(dp[int(0.5*len(t)):]), key=operator.itemgetter(1))#max(dp[int(0.5*len(t)):])
			cd_min=min(cd[int(0.5*len(t)):])
			cd_max_i, cd_max = max(enumerate(cd[int(0.5*len(t)):]), key=operator.itemgetter(1))
			cl_min=min(cl[int(0.5*len(t)):])
			cl_max_i, cl_max = max(enumerate(cl[int(0.5*len(t)):]), key=operator.itemgetter(1))

						
			#Plot lift	
			AX[count_FILE].plot(t,cl,label=self.label[count_file])
			AX[count_FILE].set_ylabel(r'$C_{l}$')
			AX[count_FILE].set_ylim(-1.1,1.1)
			AX[count_FILE].set_xlim(7.,8.)
			

			#Align figures together
			FIG[0].subplots_adjust(hspace=0)

			#Compute Period and Frequency
			cl_max_temp=cl_max
			diff_old=0
			count=0
			t0=[]
			for i in range(cl_max_i,len(cl)):
				diff=cl_max_temp-cl[i]
				cl_max_temp=cl[i]
				if diff_old<0 and diff>0:
					t0.append([i,t[i]])
					count+=1
				if count==2:break
				diff_old=diff
			period=t0[1][1]-t0[0][1]
			freq  =1./period

			#Compute pressure gradient at t0+1/2f
			index_dp=int(t0[0][0]+(t0[1][0]-t0[0][0])/2.)
			
			print "-"*100
			print self.cases[count_file]
			print "-"*100
			print  "Cd Max upper bound (Schafer & Turek) :", 3.2400
			print  "Cd Max computed                      :", cd_max
			print  "cd Max lower bound (Schafer & Turek) :", 3.2200
			print "-"*60 
			print  "Cl Max upper bound (Schafer & Turek) :", 1.0100
			print  "Cl Max computed                      :", cl_max
			print  "cl Max lower bound (Schafer & Turek) :", 0.9900
			print "-"*60 
			print  "St upper bound (Schafer & Turek) :", 0.3050
			print  "St computed                      :", (0.1*freq)/1.
			print  "St lower bound (Schafer & Turek) :", 0.2950
			print "-"*60 
			print  "dp upper bound (Schafer & Turek) :", 2.5000
			print  "dp computed                      :", dp[index_dp]
			print  "dp lower bound (Schafer & Turek) :", 2.4600

			#count_FILE+=1
			count_file+=1


CYLINDER_PLOT().plot_data()
	
plt.xlabel("Time")
plt.axhline(y=1., color='k',linestyle='--',label="Schafer & Turek")
plt.legend()
plt.show()
		
				
		






		
	


