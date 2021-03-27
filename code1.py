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
		
def Grashof_criterion(config):

	p=np.array([config['L1'], config['L2'], config['L3'], config['L4']])
	p=np.sort(p)
	sum1=p[0]+p[3]
	sum2=p[1]+p[2]
	    
	if ((sum1<sum2) and p(1)==l1):
		Type='Grashoff class-I Double Crank Mechanism'
	elif (((sum1<sum2) and p(1)==l2) | ((sum1<sum2) and p(1)==l4)):
		Type='Grashoff class-I Crank Rocker Mechanism'
	elif ((sum1<sum2) and p(1)==l3):
		Type='Grashoff class-I Double Rocker Mechanism'
	elif ((sum1<sum2) and p(1)==l4):
		Type='Grashoff class-I Rocker Crank Mechanism'
	elif (sum1>sum2):
		Type='Grashoff class-II Triple Rocker Mechanism'
	elif (sum1==sum2 and p(1)==l1 ):
		Type='Grashoff class-III Double Crank Mechanism'
	elif (sum1==sum2 and p(1)==l2 or ((sum1==sum2) and p(1)==l4 )):
		Type='Grashoff class-III Crank Rocker Mechanism'
	elif (sum1==sum2 and p(1)==l3 ):
		Type='Grashoff class-III(Deltoid Linkage) Double Rocker Mechanism'
	elif (sum1==sum2 and (l1==l2==l3==l4)):
		Type='Square with triple change point';
	elif (sum1==sum2 and (l1+l2==l3+l4) or ((sum1==sum2) and l2+l3==l1+l4 )):
		Type='Grashoff class-III Crank Rocker Mechanism'
	else:
		Type='Link lengths are incorrect.Check input'
   
	return Type
		
def points_to_FBM(config):

	config['L1'] = ((config['X3'] - config['X1'])**2 + (config['X4'] - config['X2'])**2)**0.5
	config['L2'] = ((config['X5'] - config['X3'])**2 + (config['X6'] - config['X4'])**2)**0.5
	config['L3'] = ((config['X7'] - config['X5'])**2 + (config['X8'] - config['X6'])**2)**0.5
	config['L4'] = ((config['X1'] - config['X7'])**2 + (config['X2'] - config['X8'])**2)**0.5

	config['theta1'] = np.arctan((config['X2'] - config['X4']) / (config['X1'] - config['X3']))
	config['theta2'] = np.arctan((config['X6'] - config['X4']) / (config['X5'] - config['X3']))
	config['theta3'] = np.arctan((config['X8'] - config['X6']) / (config['X7'] - config['X5']))
	config['def_theta3'] = config['theta3']
	
	config['Xa'] = (config['X1'] + config['X3'])/2
	config['Xb'] = (np.tan(config['theta1']) * (config['Xa']-config['X1'])) + config['X2']
	
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
	
	
def push_off_distance(config):
	
	angle = (config['theta3'] - config['def_theta3'])*-0.5
	cp=np.matmul(np.array([[np.cos(angle), -np.sin(angle)],
				[np.sin(angle), np.cos(angle)]]),
				np.array([config['def_Xlra'], config['def_Xlrb']]))
				
	config['Xlra'], config['Xlrb'] = cp[0], cp[1]
	
	angle = -1*angle
	
	hp=np.matmul(np.array([[np.cos(angle), -np.sin(angle)],
				[np.sin(angle), np.cos(angle)]]),
				np.array([config['def_Xhpa'], config['def_Xhpb']]))

	config['Xhpa'], config['Xhpb'] = hp[0], hp[1]

	p=((config['Xhpb']-config['Xlrb'])*config['X9']+(config['Xlra']-config['Xhpa']) \
	*config['X10']-(config['Xlra']*config['Xhpb']-config['Xlrb']*config['Xhpa']))/ \
	((config['Xhpb']-config['Xlrb'])**2+(config['Xlra']-config['Xhpa'])**2)**0.5
	
	q=((config['Xlra']-config['X9'])**2+(config['Xlrb']-config['X10'])**2-(p)**2)**0.5
	
	return p, q
		
	
def heel_contact(config):

	angle = (config['theta3'] - config['def_theta3'])*-0.5
	cp=np.matmul(np.array([[np.cos(angle), -np.sin(angle)],
				[np.sin(angle), np.cos(angle)]]),
				np.array([config['def_Xlla'], config['def_Xllb']]))
				
	config['Xlla'], config['Xllb'] = cp[0], cp[1]
	
	angle = -1*angle
	
	hp=np.matmul(np.array([[np.cos(angle), -np.sin(angle)],
				[np.sin(angle), np.cos(angle)]]),
				np.array([config['def_Xhpa'], config['def_Xhpb']]))

	config['Xhpa'], config['Xhpb'] = hp[0], hp[1]

	p=((config['Xhpb']-config['Xllb'])*config['X9']+(config['Xlla']-config['Xhpa']) \
	*config['X10']-(config['Xlla']*config['Xhpb']-config['Xllb']*config['Xhpa']))/ \
	((config['Xhpb']-config['Xllb'])**2+(config['Xlla']-config['Xhpa'])**2)**0.5
	
	q=((config['Xlla']-config['X9'])**2+(config['Xllb']-config['X10'])**2-(p)**2)**0.5
	
	return p, q
		
	
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


def plot(config, i):

	l1=ax[0].plot([config['X1'], config['X3']],[config['X2'], config['X4']],linestyle='-', marker='x', color='b',label='l1')
	l2=ax[0].plot([config['X3'], config['X5']],[config['X4'], config['X6']],linestyle='-', marker='x', color='g',label='l2')
	l3=ax[0].plot([config['X5'], config['X7']],[config['X6'], config['X8']],linestyle='-', marker='x', color='r',label='l3')
	l4=ax[0].plot([config['X7'], config['X1']],[config['X8'], config['X2']],linestyle='-', marker='x', color='c',label='l4')
	l5=ax[0].plot([config['X5'], config['X9']],[config['X6'], config['X10']],linestyle='--', color='k')
	l6=ax[0].plot([config['X7'], config['X9']],[config['X8'], config['X10']],linestyle='--', color='k')
	femur = ax[0].plot([config['Xua'], config['X7']],[config['Xub'], config['X8']],linestyle='-', color='b',label='femur')
	ll1 = ax[0].plot([config['Xlla'], config['Xua']],[config['Xllb'], config['Xub']],linestyle='--', color='r',label='anterior load line')
	ll2 = ax[0].plot([config['Xlra'], config['Xua']],[config['Xlrb'], config['Xub']],linestyle='--', color='r',label='anterior load line')
	ax[0].legend(loc='best')
	cp=ax[0].scatter(config['X9'], config['X10'], marker='.', color='k')
	p1, q1 = push_off_distance(config)
	p2, q2 = heel_contact(config)
	ax[1].scatter(i, p2/q2, marker='.', color='b')
	ax[1].scatter(i, p1/q1, marker='.', color='k')
	ax[1].plot([0, i+1], [0, 0], color='k')
	plt.pause(config['pause'])
	l1.pop(0).remove()
	l2.pop(0).remove()
	l3.pop(0).remove()
	l4.pop(0).remove()
	l5.pop(0).remove()
	l6.pop(0).remove()
	femur.pop(0).remove()
	ll1.pop(0).remove()
	ll2.pop(0).remove()

config = {
		'X1':np.nan,'X2':np.nan,'X3':np.nan,'X4':np.nan,
		'X5':np.nan,'X6':np.nan,'X7':np.nan,'X8':np.nan,
		'Xa':np.nan, 'Xb': np.nan, 'Xlma':0, 'Xlmb':-500,
		'Xlla':-50, 'Xllb':-500,'Xlra':200, 'Xlrb':-500,
		'Xua':0, 'Xub':400,'Xe':0, 'Xf':400,
		'def_Xlla':-50, 'def_Xllb':-500,
		'def_Xlra':200, 'def_Xlrb':-500,
		'Xhpa':0, 'Xhpb':400,
		'def_Xhpa':0, 'def_Xhpb':400,
		'L1':np.nan,'L2':np.nan,'L3':np.nan,'L4':np.nan,
		'theta1':np.nan, 'theta2':np.nan,
		'theta3':np.nan, 'theta4':np.nan,
		'def_theta3':np.nan, 'stance':True,
		'draw':1,'fc': 10, 'pause':.2, 'inc_fac': +180,
		'frames':40, 'filename':"points.csv"
		}
		

read_file(config)
points_to_FBM(config)
Type = Grashof_criterion(config)
print(colored(Type, 'blue'))
fig,ax = plt.subplots(1,2, figsize=(15,10))
feet = ax[0].plot([config['Xlla'], config['Xlra']],[config['Xllb'], config['Xlrb']],linestyle='-', color='b',label='feet')
Tibia = ax[0].plot([config['Xlma'], config['Xa']],[config['Xlmb'], config['Xb']],linestyle='-', color='b',label='Tibia')
mx = max(config['L1'], config['L2'], config['L3'], config['L4'])
fc = config['fc']
ax[0].set_xlim(-fc*mx, fc*mx)
ax[0].set_ylim(-fc*mx, fc*mx)
ax[1].scatter(None, None, marker='.', color='b', label='Heel Contact')
ax[1].scatter(None, None, marker='.', color='k', label='Push Off')
ax[1].legend(loc='best')

for i in range(config['frames']):
	_FBM_angles(config)
	if not config['draw']:
		print(colored(f"Error at theta3 = %f" % config['theta3'], 'red'))
		config['theta3'] += np.pi/config['inc_fac']
		continue
	FBM_simulations(config)
	_FBM_centroid(config)
	plot(config, i)
	#print("Center = ",(config['X9'], config['X10']))
	#print("Knee Center = ", (config['X7'], config['X8']), "\n")
	config['theta3'] += np.pi/config['inc_fac']

