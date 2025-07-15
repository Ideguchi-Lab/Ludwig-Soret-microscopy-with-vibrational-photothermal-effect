
# %%
import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
import os
import tqdm

path = '/home/toda/soret/20220903_cell7_nuc/'
os.chdir(path)

# parameter setting diffusion calculation
DX = 0.207 # pixel size (um)
DY = 0.207
NX = 141
NY = 141
DT = 1000 # step of time (us)
NT = int(2.4*10**3 + 1) # the number of step
T_OUT = 296 # boundary temperature
N_SAVE = 100 # saving step interval

Tij = cp.asarray(cp.loadtxt("temperature.csv", delimiter=',', dtype=cp.float64))
mask = cp.asarray(cp.loadtxt("mask_crop.csv", delimiter=',', dtype=cp.float64))

mask_edge_x = cp.ones((NX, NY))
mask_edge_x[1:-1,1:-1] = (mask[2:,1:-1]-mask[1:-1,1:-1])
#mask_edge_x[1:-1,1:-1] *= ((mask_edge_x[1:-1,1:-1] > 0) + (mask_edge_x [2:,1:-1]) < 0)
mask_edge_y = cp.ones((NX, NY))
mask_edge_y[1:-1,1:-1] = (mask[1:-1,2:]-mask[1:-1,1:-1])
#mask_edge_y[1:-1,1:-1] *= ((mask_edge_y[1:-1,1:-1] > 0) + (mask_edge_y[1:-1,2:]) < 0) 

C_n = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))
C_0 = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))

Tij = Tij * (mask > 0)
Tij2 = cp.zeros((NX, NY))

C_02 = cp.zeros((NX, NY))
C_n2 = cp.zeros((NX, NY))

#print(cp.sum(C_n*mask,[0,1]))

Tij2[1:-1,1:-1] = (
        (mask_edge_x[1:-1,1:-1] == 1)*Tij[2:,1:-1] 
        + (mask_edge_x[:-2,1:-1] == -1)*Tij[:-2,1:-1] 
        + (1-(mask_edge_x[1:-1,1:-1] == 1))* (1-(mask_edge_x[:-2,1:-1] == -1))*(mask_edge_y[1:-1,1:-1] == 1)*Tij[1:-1,2:] 
        + (1-(mask_edge_x[1:-1,1:-1] == 1))* (1-(mask_edge_x[:-2,1:-1] == -1))*(mask_edge_y[1:-1,:-2] == -1)*Tij[1:-1,:-2]
        + (mask[1:-1,1:-1] > 0)*Tij[1:-1,1:-1]* (1-(mask_edge_x[1:-1,1:-1] == 1))*(1-(mask_edge_y[1:-1,1:-1] == 1))*(1-(mask_edge_x[:-2,1:-1] == -1))*(1-(mask_edge_y[1:-1,:-2] == -1))
    )

C_n2[1:-1,1:-1] = (
        (mask_edge_x[1:-1,1:-1] == 1)*C_n[2:,1:-1] 
        + (mask_edge_x[:-2,1:-1] == -1)*C_n[:-2,1:-1] 
        + (1-(mask_edge_x[1:-1,1:-1] == 1))* (1-(mask_edge_x[:-2,1:-1] == -1))*(mask_edge_y[1:-1,1:-1] == 1)*C_n[1:-1,2:] 
        + (1-(mask_edge_x[1:-1,1:-1] == 1))* (1-(mask_edge_x[:-2,1:-1] == -1))*(mask_edge_y[1:-1,:-2] == -1)*C_n[1:-1,:-2]
        + (mask[1:-1,1:-1] > 0)*C_n[1:-1,1:-1]* (1-(mask_edge_x[:-2,1:-1] == -1))*(1-(mask_edge_y[1:-1,:-2] == -1))
        + (mask[1:-1,1:-1] == 0)*1* (1-(mask_edge_x[1:-1,1:-1] == 1))*(1-(mask_edge_y[1:-1,1:-1] == 1))*(1-(mask_edge_x[:-2,1:-1] == -1))*(1-(mask_edge_y[1:-1,:-2] == -1))
    )

C_02[1:-1,1:-1] = (
        (mask_edge_x[1:-1,1:-1] == 1)*C_0[2:,1:-1] 
        + (mask_edge_x[:-2,1:-1] == -1)*C_0[:-2,1:-1] 
        + (1-(mask_edge_x[1:-1,1:-1] == 1))* (1-(mask_edge_x[:-2,1:-1] == -1))*(mask_edge_y[1:-1,1:-1] == 1)*C_0[1:-1,2:] 
        + (1-(mask_edge_x[1:-1,1:-1] == 1))* (1-(mask_edge_x[:-2,1:-1] == -1))*(mask_edge_y[1:-1,:-2] == -1)*C_0[1:-1,:-2]
        + (mask[1:-1,1:-1] > 0)*C_0[1:-1,1:-1]* (1-(mask_edge_x[:-2,1:-1] == -1))*(1-(mask_edge_y[1:-1,:-2] == -1))
        + (mask[1:-1,1:-1] == 0)*1* (1-(mask_edge_x[1:-1,1:-1] == 1))*(1-(mask_edge_y[1:-1,1:-1] == 1))*(1-(mask_edge_x[:-2,1:-1] == -1))*(1-(mask_edge_y[1:-1,:-2] == -1))
    )

# %%

def ftcs(Cij,dt1,d1):
    Cij_copy = Cij.copy()
    Cij_copy[2:-2,2:-2] += dt1*DT*Cij_copy[2:-2,2:-2]*(
        Tij2[3:-1,2:-2] + Tij2[1:-3,2:-2]
        +Tij2[2:-2,3:-1] + Tij2[2:-2,1:-3]
        -4*Tij2[2:-2,2:-2])/(DX*DY) + d1*DT*(
        (Cij_copy[3:-1,2:-2]-C_02[3:-1,2:-2]) + (Cij_copy[1:-3,2:-2]-C_02[1:-3,2:-2])
        + (Cij_copy[2:-2,3:-1]-C_02[2:-2,3:-1]) + (Cij_copy[2:-2,1:-3]-C_02[2:-2,1:-3]) 
        -4*(Cij_copy[2:-2,2:-2]-C_02[2:-2,2:-2]))/(DX*DY) + dt1*DT*(
        (Tij2[3:-1,2:-2]-Tij2[2:-2,2:-2])/DX * (Cij_copy[3:-1,2:-2]-Cij_copy[2:-2,2:-2])/DX
        + (Tij2[2:-2,3:-1]-Tij2[2:-2,2:-2])/DY * (Cij_copy[2:-2,3:-1]-Cij_copy[2:-2,2:-2])/DY
        )  

    Cij_copy[:,:] = Cij_copy[:,:]*(mask > 0) 
    
    Cij_copy[2:-2,2:-2] = ((mask_edge_x[1:-3,2:-2] == 1)*(Cij_copy[3:-1,2:-2]-(C_02[3:-1,2:-2]-C_02[2:-2,2:-2]))/(1 - dt1/d1*(Tij2[3:-1,2:-2]-Tij2[2:-2,2:-2]))
        + (mask_edge_x[2:-2,2:-2] == -1)*(Cij_copy[1:-3,2:-2]*(1 - dt1/d1*(Tij2[2:-2,2:-2]-Tij2[1:-3,2:-2])) + (C_02[2:-2,2:-2]-C_02[1:-3,2:-2])) 
        + (1-(mask_edge_x[1:-3,2:-2] == 1))* (1-(mask_edge_x[2:-2,2:-2] == -1))*(mask_edge_y[2:-2,1:-3] == 1)*(Cij_copy[2:-2,3:-1]-(C_02[2:-2,3:-1]-C_02[2:-2,2:-2]))/(1 - dt1/d1*(Tij2[2:-2,3:-1]-Tij2[2:-2,2:-2])) 
        + (1-(mask_edge_x[1:-3,2:-2] == 1))* (1-(mask_edge_x[2:-2,2:-2] == -1))*(mask_edge_y[2:-2,2:-2] == -1)*(Cij_copy[2:-2,1:-3]*(1 - dt1/d1*(Tij2[2:-2,2:-2]-Tij2[2:-2,1:-3])) +  (C_02[2:-2,2:-2]-C_02[2:-2, 1:-3])) 
        + (mask[2:-2,2:-2] > 0)*(1-(mask_edge_x[2:-2,2:-2] == -1))*(1-(mask_edge_y[2:-2,2:-2] == -1))* (1-(mask_edge_x[1:-3,2:-2] == 1))*(1-(mask_edge_y[2:-2,1:-3] == 1))*Cij_copy[2:-2,2:-2]
    )

    Cij_copy[2:-2,2:-2] += (mask[2:-2,2:-2] == 0)*(
        (mask_edge_x[2:-2,2:-2] == 1)*Cij_copy[3:-1,2:-2]  
        + (mask_edge_x[1:-3,2:-2] == -1)* Cij_copy[1:-3,2:-2]   
        + (1-(mask_edge_x[2:-2,2:-2] == 1))* (1-(mask_edge_x[1:-3,2:-2] == -1))* (mask_edge_y[2:-2,2:-2] == 1) *Cij_copy[2:-2,3:-1]
        + (1-(mask_edge_x[2:-2,2:-2] == 1))* (1-(mask_edge_x[1:-3,2:-2] == -1))* (mask_edge_y[2:-2,1:-3] == -1)* Cij_copy[2:-2,1:-3]
        + (mask[2:-2,2:-2] == 0)* (1-(mask_edge_x[2:-2,2:-2] == 1))* (1-(mask_edge_x[1:-3,2:-2] == -1))* (1-(mask_edge_y[2:-2,2:-2] == 1))*(1-(mask_edge_y[2:-2,1:-3] == -1))*1
    )
   
    return Cij_copy


obj = cp.ones((NX, NY, NT//N_SAVE+1 )) 

D1 = 4.4945*10**-6 # diffusion coefficient [um^2/us]
St = 0.0092874 #0.0092 #20220903 cell7 nuc
Dt1 = D1*St # thermophoretic mobility [um^2/us]
xcenter=70
ycenter=70
temporal_evolution = cp.asarray(cp.loadtxt("temporal_evolution_drymass.csv", delimiter=',', dtype=cp.float64))
time= np.arange(0, (NT-1)//N_SAVE+1)*N_SAVE*DT*10**-6
a=temporal_evolution[5:127:5,1]
b=temporal_evolution[5:127:5,0]
#print(a, a.shape)
#print(b, b.shape)

C_0 = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))

def ff(D1, St):
    Dt1 = D1 * St  # thermophoretic mobility
    C_n = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))

    for n in range(NT):
        C_np1 = cp.ones((NY, NX)) 
        C_np1 = ftcs(C_n, Dt1, D1)

        if (n % N_SAVE == 0):
            #print(n)
            obj[:, :, n // N_SAVE] = C_n[:, :]
        C_n = C_np1

    temp = cp.mean(obj[xcenter-5:xcenter+5, ycenter-5:ycenter+5, :], axis=(0, 1))-cp.mean(C_0[xcenter-5:xcenter+5, ycenter-5:ycenter+5])
    return temp


# %%
def fitting_function(params):
    D1, St = params
    temp = ff(D1, St)
    # 目標のデータとの差の二乗和を返す
    return cp.sum(cp.abs(temp+a) ** 2).item()

from scipy.optimize import minimize

# 初期値を設定
initial_guess = [4.6e-6, 0.0093]  # D1, St の初期値

# 最小化を実行
result = minimize(fitting_function, initial_guess, method='Nelder-Mead')

# 結果を表示
fitted_D1, fitted_St = result.x
print("Fitted D1:", fitted_D1)
print("Fitted St:", fitted_St)
# print("Residual Sum of Squares:", result.fun)

# 例: 初期条件ファイルのロード
C_0 = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))

def fit_ffff(x, D1, St):
    return ffff(x, D1, St)  # `ffff`関数自体の戻り値は `np.array`


def ffff(x, D1, St):
    Dt1 = D1 * St  # thermophoretic mobility
    C_n = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))
    
    # 保存するタイミングを決定
    max_n = int(x.max() * N_SAVE)
    print(max_n)
    obj = cp.zeros((NY, NX, x.size))  # 必要なxの数だけ保存用配列を用意

    index = 0  # 保存用インデックス
    for n in tqdm(range(max_n)):
        C_np1 = ftcs(C_n, Dt1, D1)  # 計算（関数は適宜変更）
        
        # 指定した保存タイミングでのみ結果を保存
        if n % N_SAVE == 0 and (n // N_SAVE) in x:
            obj[:, :, index] = C_n[:, :]
            index += 1
        
        C_n = C_np1  # 状態更新

    # 中心領域の平均を計算し、tempを返す
    temp = cp.mean(obj[xcenter-5:xcenter+5, ycenter-5:ycenter+5, :], axis=(0, 1)) - cp.mean(C_0[xcenter-5:xcenter+5, ycenter-5:ycenter+5])
    return temp.get()  # cupy配列からnumpy配列に変換して返す


# 誤差の範囲を計算するために、最適化前の状態から少しずつ変化させてみる
delta_d1 = 0.02e-6  # 誤差幅を決定する値
delta_st = 0.02e-3  # 誤差幅を決定する値
D1_errors = []
St_errors = []

for d1_variation in np.arange(-delta_d1*5, delta_d1*5, delta_d1):
    for st_variation in np.arange(-delta_st*5, delta_st*5, delta_st):
        test_D1 = fitted_D1 + d1_variation
        test_St = fitted_St + st_variation
        temp_test = ff(test_D1, test_St)
        test_dev = cp.sum(cp.abs(temp_test + a)** 2) 
        
        if test_dev <= result.fun * 1.1:  # 10%増加に収まる場合
            D1_errors.append(test_D1)
            St_errors.append(test_St)

print (D1_errors)
print("D1 Error Range:", min(D1_errors), max(D1_errors))
print("St Error Range:", min(St_errors), max(St_errors))

# %%
def fff(D1, St):
    Dt1 = D1 * St  # thermophoretic mobility
    C_n = cp.asarray(cp.loadtxt("initial_con_drymass.csv", delimiter=',', dtype=cp.float64))

    for n in range(NT):
        C_np1 = cp.ones((NY, NX)) 
        C_np1 = ftcs(C_n, Dt1, D1)

        if (n % N_SAVE == 0):
            obj[:, :, n // N_SAVE] = (C_n[:, :]-C_0)*mask
        C_n = 1 * C_np1

    return obj

temp = ff(fitted_D1, fitted_St)
sim = fff(fitted_D1, fitted_St)
sim *= mask[:,:,np.newaxis]

plt.plot(time, cp.asnumpy(-a), 'o')
plt.plot(time, cp.asnumpy(temp),'o')
#plt.plot(time, cp.asnumpy(glass))
plt.show()

plt.imshow((cp.asnumpy(sim[:,:,24])))
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar()
plt.show()

#%%
for t_scan in range(0, 25):
    file_name = f'/home/toda/soret/20220903_cell7_nuc/sim_{t_scan}s.csv'
    cp.savetxt(file_name, sim[:,:,t_scan])    

sim_temporal = temp
file_name = f'/home/toda/soret/20220903_cell7_nuc/sim_temporal.csv'
cp.savetxt(file_name, temp)

# %%
