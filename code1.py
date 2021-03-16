import numpy as np
from termcolor import colored
import matplotlib.pyplot as plt
import pandas as pd

def knee_center():

	c1 = 0.0791
	c2 = 5.733E-4
	c3 = 7.682E-6
	c4 = 5.759E-8

	c5 = 0.3695
	c6 = 2.958E-3
	c7 = 7.666E-6

	c8 = -0.0683
	c9 = 8.804E-4
	c10 = 3.750E-6

	c11 = -0.1283
	c12 = 4.796E-4

	for Beta in range(0, 120):

		Varus = (c1 * Beta) - (c2 * Beta**2) - (c3 * Beta**3) + (c4 * Beta**4)
		IntRot = (c5 * Beta) - (c6 * Beta**2) + (c7 * Beta**3)
		ydist = (c8 * Beta) + (c9 * Beta**2) - (c10 * Beta**3)
		zdist = (c11 * Beta) + (c12 * Beta**2)

		X1 = -60
		Y1 = -np.sin(Varus) * X1 + ydist
		Z1 = np.cos(Varus) * np.sin(IntRot) * X1 + zdist
		
def read_file(config):

	FBM_type = int(input("Enter type number 0 - 6 \n"))
	df = pd.read_csv(config['filename'], header = None)
	co_ordinates = 	df.iloc[FBM_type,:]
	co_ordinates = [float(point) for point in co_ordinates]
	for number, point in enumerate(co_ordinates):
		config[f'X%d'%(number+1)] = point
		
def points_to_FBM(config):

	config['L1'] = ((config['X3'] - config['X1'])**2 + (config['X4'] - config['X2'])**2)**0.5
	config['L2'] = ((config['X5'] - config['X3'])**2 + (config['X6'] - config['X4'])**2)**0.5
	config['L3'] = ((config['X7'] - config['X5'])**2 + (config['X8'] - config['X6'])**2)**0.5
	config['L4'] = ((config['X1'] - config['X7'])**2 + (config['X2'] - config['X8'])**2)**0.5

	config['theta1'] = np.arctan((config['X2'] - config['X4']) / (config['X1'] - config['X3']))
	config['theta2'] = np.arctan((config['X6'] - config['X4']) / (config['X5'] - config['X3']))
	config['theta3'] = np.arctan((config['X8'] - config['X6']) / (config['X7'] - config['X5']))
	
	print(f"L1=%0.2f, L2=%0.2f, L3=%0.2f, L4=%0.2f \ntheta1=%0.2f, theta2=%0.2f, theta3=%0.2f"%
	(config['L1'], config['L2'], config['L3'], config['L4'], 180/np.pi * config['theta1'], 
	180/np.pi * config['theta2'], 180/np.pi * config['theta3']))

	
def _FBM_angles(config):

	L1, L2, L3, L4, theta1, theta3 = config['L1'],config['L2'],\
	 config['L3'],config['L4'],config['theta1'], config['theta3']

	K = (L1**2 + L2**2 + L3**2 - L4**2)/2
	A = K - L2*L3*np.cos(theta3) + L1*L2*np.cos(theta1) - L1*L3*np.cos(theta1-theta3)
	B = 2*(L2*L3*np.sin(theta3) - L1*L2*np.sin(theta1))
	C = K + L2*L3*np.cos(theta3) - L1*L2*np.cos(theta1) - L1*L3*np.cos(theta1-theta3)
	
	if (B**2 - 4*A*C)>=0:
		theta21 = 2*np.arctan((-B+(B**2 - 4*A*C)**0.5)/(2*A))
		theta22 = 2*np.arctan((-B-(B**2 - 4*A*C)**0.5)/(2*A))
	
	else:
		config['draw']=0
		return config
	
	if -1<=((L2*np.sin(theta21)+ L3*np.sin(theta3) - L1*np.sin(theta1))/L4)<=1:
		theta4 = np.arcsin((L2*np.sin(theta21)+ L3*np.sin(theta3) - L1*np.sin(theta1))/L4)
		theta2 = theta21
	else:
		theta4 = np.arcsin((L2*np.sin(theta22)+ L3*np.sin(theta3) - L1*np.sin(theta1))/L4)
		theta2 = theta22
		
	config['theta2']=theta2
	config['theta4']=theta4
		
	
def _FBM_centroid(config):

	A = np.linalg.inv(np.array([[(config['X6']-config['X4']),(config['X3']-config['X5'])],
					[(config['X8'] - config['X2']),(config['X1']-config['X7'])]]))
	B = np.array([config['X3']*config['X6'] - config['X4']*config['X5'],
					config['X1']*config['X8'] - config['X2']*config['X7']])
	C = np.matmul(A,B)
	
	config['X9'] = C[0]
	config['X10'] = C[1]
	
def FBM_simulations(config):


	if config['stance']:
	
		config['X5'] =  config['X3'] + config['L2']*np.cos(config['theta2'])
		config['X6'] =  config['X4'] + config['L2']*np.sin(config['theta2'])
		
		config['X7'] =  config['X5'] + config['L3']*np.cos(config['theta3'])
		config['X8'] =  config['X6'] + config['L3']*np.sin(config['theta3'])

	else:

		config['X1'] =  config['X7'] - config['L4']*np.cos(config['theta4'])
		config['X2'] =  config['X8'] - config['L4']*np.sin(config['theta4'])

		config['X3'] =  config['X1'] - config['L1']*np.cos(config['theta1'])
		config['X4'] =  config['X2'] - config['L1']*np.sin(config['theta1'])


def plot(config):

	l1=ax.plot([config['X1'], config['X3']],[config['X2'], config['X4']],linestyle='-', marker='x', color='b',label='l1')
	l2=ax.plot([config['X3'], config['X5']],[config['X4'], config['X6']],linestyle='-', marker='x', color='g',label='l2')
	l3=ax.plot([config['X5'], config['X7']],[config['X6'], config['X8']],linestyle='-', marker='x', color='r',label='l3')
	l4=ax.plot([config['X7'], config['X1']],[config['X8'], config['X2']],linestyle='-', marker='x', color='c',label='l4')
	l5=ax.plot([config['X5'], config['X9']],[config['X6'], config['X10']],linestyle='--', color='k')
	l6=ax.plot([config['X7'], config['X9']],[config['X8'], config['X10']],linestyle='--', color='k')
	ax.legend(loc='best')
	cp=ax.scatter(config['X9'], config['X10'], marker='.', color='k')
	plt.pause(config['pause'])
	l1.pop(0).remove()
	l2.pop(0).remove()
	l3.pop(0).remove()
	l4.pop(0).remove()
	l5.pop(0).remove()
	l6.pop(0).remove()
	

config = {
		'X1':np.nan,'X2':np.nan,'X3':np.nan,'X4':np.nan,
		'X5':np.nan,'X6':np.nan,'X7':np.nan,'X8':np.nan,
		'L1':np.nan,'L2':np.nan,'L3':np.nan,'L4':np.nan,
		'theta1':np.nan, 'theta2':np.nan,
		'theta3':np.nan, 'theta4':np.nan,
		'stance':True,
		'draw':1,'fc': 3, 'pause':0.01, 'inc_fac': +180,
		'frames':120, 'filename':"points.csv"
		}
		

read_file(config)
points_to_FBM(config)
fig,ax = plt.subplots(1,1, figsize=(10,10))
mx = max(config['L1'], config['L2'], config['L3'], config['L4'])
fc = config['fc']
plt.xlim(-fc*mx, fc*mx)
plt.ylim(-fc*mx, fc*mx)

for i in range(config['frames']):
	_FBM_angles(config)
	if not config['draw']:
		print(colored(f"Error at theta3 = %f" % config['theta3'], 'red'))
		config['theta3'] += np.pi/config['inc_fac']
		continue
	FBM_simulations(config)
	_FBM_centroid(config)
	plot(config)
	print("Center = ",(config['X9'], config['X10']))
	config['theta3'] += np.pi/config['inc_fac']

