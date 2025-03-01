import numpy as np
import matplotlib.pyplot as plt
from random import shuffle
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
import tkinter as tk
from tkinter import ttk

def monte_carlo(mu, sigma, k=0, IC = 999, nb=100, taille=1000):
    L_UB = []
    L_LB = []
    if taille < 1000:
        taille = 1000
    for w in range(nb):
        L_donnees = np.random.normal(loc = mu, scale=3*sigma, size=taille)
        L_C_p = [0]
        L_C_m = [0]
        for i in range(len(L_donnees)):
            L_C_p.append(max(0, L_C_p[i] + L_donnees[i] - (mu - k)))
            L_C_m.append(min(0, L_C_m[i] + L_donnees[i] - (mu - k)))
        L_C_p.sort()
        L_C_m.sort()

        L_LB.append(L_C_m[IC])
        L_UB.append(L_C_p[1000-IC])
    return np.mean(L_UB), np.mean(L_LB)

def cusum_p(L_donnees, k=0):
    """
    CONTINUOUS INSPECTION SCHEMES Get access Arrow, E. S. Page, Biometrika, 
    Volume 41, Issue 1-2, June 1954, Pages 100–115, https://doi.org/10.1093/biomet/41.1-2.100
    """
    moy = np.mean(L_donnees)
    L_C = [0]
    for i in range(len(L_donnees)):
        L_C.append(max(0, L_C[i] + L_donnees[i] - (moy - k)))
    return L_C

def bootstrap_cusum_p(L_donnees, k=0, M = 1000, IC = 0.95, estimateur=0):
    L_C = cusum_p(L_donnees)
    moy = np.mean(L_donnees)
    C_diff = np.max(L_C)-np.min(L_C)
    compt = 0
    L_donnees_bootstrap = [elmt for elmt in L_donnees]
    # Bootstrap :
    for i in range(M):
        shuffle(L_donnees_bootstrap)
        L_C_bootstrap = [0]
        for i in range(len(L_donnees_bootstrap)):
            L_C_bootstrap.append(max(0, L_C_bootstrap[i] + L_donnees_bootstrap[i] - (moy - k)))

        if np.max(L_C_bootstrap)-np.min(L_C_bootstrap) < C_diff:
            compt += 1
    ind = -1
    if compt/M >= IC/100 :
        MSE_m = []
        if len(L_donnees)> 3 and estimateur==1:
            for i in range(1, len(L_donnees)-2):
                MSE_m.append(np.sum(np.power(L_donnees[:i]-np.mean(L_donnees[:i]), 2))\
                             +np.sum(np.power(L_donnees[i:]-np.mean(L_donnees[i:]), 2)))
            ind = np.argmin(MSE_m)+1
        else:
            ind = np.argmax(L_C)
    return ind

def cusum_m(L_donnees, k=0):
    """
    CONTINUOUS INSPECTION SCHEMES Get access Arrow, E. S. Page, Biometrika, 
    Volume 41, Issue 1-2, June 1954, Pages 100–115, https://doi.org/10.1093/biomet/41.1-2.100
    """
    moy = np.mean(L_donnees)
    L_C = [0]
    for i in range(len(L_donnees)):
        L_C.append(min(0, L_C[i] + L_donnees[i] - (moy - k)))
    return L_C

def bootstrap_cusum_m(L_donnees, k=0, M = 1000, IC = 0.95, estimateur=0):
    L_C = cusum_m(L_donnees)
    moy = np.mean(L_donnees)
    C_diff = np.max(L_C)-np.min(L_C)
    compt = 0
    L_donnees_bootstrap = [elmt for elmt in L_donnees]
    # Bootstrap :
    for i in range(M):
        shuffle(L_donnees_bootstrap)
        L_C_bootstrap = [0]
        for i in range(len(L_donnees_bootstrap)):
            L_C_bootstrap.append(min(0, L_C_bootstrap[i] + L_donnees_bootstrap[i] - (moy - k)))

        if np.max(L_C_bootstrap)-np.min(L_C_bootstrap) < C_diff:
            compt += 1
    ind = -1
    if compt/M >= IC/100 :
        MSE_m = []
        if len(L_donnees)> 3 and estimateur==1:
            for i in range(1, len(L_donnees)-2):
                MSE_m.append(np.sum(np.power(L_donnees[:i]-np.mean(L_donnees[:i]), 2))\
                             +np.sum(np.power(L_donnees[i:]-np.mean(L_donnees[i:]), 2)))
            ind = np.argmin(MSE_m)+1
        else:
            ind = np.argmin(L_C)
    return ind

class app:
    def __init__(self):
        # Variable
        self._L_ind = []
        self._L_mesures = []

        self._fn = tk.Tk()
        self._fn.title("CUSUM")
        for i in range(3):
            self._fn.rowconfigure(i, weight=1)
            self._fn.columnconfigure(i, weight=1)

        #Canvas display
        self._cadre_data()
        self._cadre_hyper_parametre()
        self._cadre_regle()
        bout = ttk.Button(self._fn, text="Calculation", width=3, 
                          command=self._calcul)
        bout.grid(row=3, column=0, columnspan=2, sticky="NEWS")

        self._s = ttk.Style()
        self._s.theme_use('winnative')

        self._fn.mainloop()

    def _cadre_data(self):
        cdr = tk.LabelFrame(self._fn, text="Data:")
        cdr.grid(row=0, column=0)
        
        self._var_nb_mesures = tk.Variable(self._fn, value =0)
        lb = ttk.Label(cdr, text = "Measures:")
        lb.grid(row=0, column=0 )
        sp = ttk.Entry(cdr, width=20,
                         textvariable = self._var_nb_mesures)
        sp.grid(row=0, column=1)

        bout = ttk.Button(cdr, text="Add",
                          command=self._affichage_variable)
        bout.grid(row=0, column=2, columnspan=2, sticky="NEWS")

        self._list_mesures = tk.Listbox(cdr, width=40, height=11,
                                         selectmode = tk.SINGLE)
        self._list_mesures.grid(row=1, column=0, columnspan=2, rowspan=3)

        barre_deroulante = ttk.Scrollbar(cdr, orient=tk.VERTICAL, command=self._list_mesures.yview)
        barre_deroulante.grid(row=1, column=2, rowspan=3, sticky='ns')

        self._list_mesures.config(yscrollcommand = barre_deroulante.set)

        # SUPP
        bout = ttk.Button(cdr, text="✕", width=3, 
                          command=self._sup)
        bout.grid(row=2, column=3, sticky="NEWS")

        # UP
        bout = ttk.Button(cdr, text="↑",
                          command=self._haut)
        bout.grid(row=1, column=3, sticky="NEWS")

        # DOWN
        bout = ttk.Button(cdr, text="↓",
                          command=self._bas)
        bout.grid(row=3, column=3, sticky="NEWS")

    def _affichage_variable(self):
        self._list_mesures.insert(tk.END, float(self._var_nb_mesures.get()))

    def _haut(self):
        for i in self._list_mesures.curselection():
            if i > 0:
                item = self._list_mesures.get(i)
                self._list_mesures.delete(i)
                self._list_mesures.insert(i-1, item)
    
    def _sup(self):
        for i in self._list_mesures.curselection():
            self._list_mesures.delete(i)

    def _bas(self):
        for i in self._list_mesures.curselection():
            if i < self._list_mesures.index("end"):
                item = self._list_mesures.get(i)
                self._list_mesures.delete(i)
                self._list_mesures.insert(i+1, item)

    def _cadre_hyper_parametre(self):
        cdr = tk.LabelFrame(self._fn, text="Hyperparameters:")
        cdr.grid(row=0, column=1)

        cdr_bootstrap = tk.LabelFrame(cdr, text="Bootstrap:")
        cdr_bootstrap.grid(row=0, column=0)

        self._var_IC = tk.Variable(self._fn, value =95)
        lb = ttk.Label(cdr_bootstrap, text = "Confidence interval (%): ")
        lb.grid(row=0, column=0)
        sp = ttk.Entry(cdr_bootstrap, width=5,
                         textvariable = self._var_IC)
        sp.grid(row=0, column=1)

        self._var_nb_bootstrap = tk.Variable(self._fn, value =1000)
        lb = ttk.Label(cdr_bootstrap, text = "Number of analyses: ")
        lb.grid(row=1, column=0, sticky='W')
        sp = ttk.Entry(cdr_bootstrap, width=5,
                         textvariable = self._var_nb_bootstrap)
        sp.grid(row=1, column=1)

        self._var_esti_bootstrap = tk.IntVar(self._fn, 1)
        lb = ttk.Label(cdr_bootstrap, text = "Estimator : ")
        lb.grid(row=2, column=0, sticky='W')
        rad = tk.Radiobutton(cdr_bootstrap, 
               text="S", variable=self._var_esti_bootstrap, 
               value=0)
        rad.grid(row=2, column=1, sticky='W')
        rad = tk.Radiobutton(cdr_bootstrap, 
               text="MSE", variable=self._var_esti_bootstrap, 
               value=1)
        rad.grid(row=3, column=1, sticky='W')

        cdr_monte_carlo = tk.LabelFrame(cdr, text="Monte-Carlo:")
        cdr_monte_carlo.grid(row=1, column=0)

        self._var_monte_carlo_IC = tk.Variable(self._fn, value =90)
        lb = ttk.Label(cdr_monte_carlo, text = "Confidence interval (%): ")
        lb.grid(row=0, column=0)
        sp = ttk.Entry(cdr_monte_carlo, width=5,
                         textvariable = self._var_monte_carlo_IC)
        sp.grid(row=0, column=1)

        self._var_nb_monte_carlo = tk.Variable(self._fn, value =100)
        lb = ttk.Label(cdr_monte_carlo, text = "Number of analyses: ")
        lb.grid(row=1, column=0, sticky='W')
        sp = ttk.Entry(cdr_monte_carlo, width=5,
                         textvariable = self._var_nb_monte_carlo)
        sp.grid(row=1, column=1)

        self._var_nb_taille_monte_carlo = tk.Variable(self._fn, value =1000)
        lb = ttk.Label(cdr_monte_carlo, text = "Sample size: ")
        lb.grid(row=2, column=0, sticky='W')
        sp = ttk.Entry(cdr_monte_carlo, width=5,
                         textvariable = self._var_nb_taille_monte_carlo)
        sp.grid(row=2, column=1)

    def _cadre_regle(self):
        cdr = tk.Frame(self._fn)
        cdr.grid(row=1, column=0, columnspan=2)

        cd = tk.LabelFrame(cdr, text="Rules:")
        cd.grid(row=0, column=0)

        self._var_regle_1 = tk.IntVar(self._fn, value=1)
        self._var_regle_2 = tk.IntVar(self._fn, value=1)
        self._var_regle_3 = tk.IntVar(self._fn, value=1)
        self._var_regle_4 = tk.IntVar(self._fn, value=1)
        C = tk.Checkbutton(cd, text = "Rule 1", variable = self._var_regle_1, \
                                onvalue = 1, offvalue = 0)
        C.grid(row=0, column=0)
        C = tk.Checkbutton(cd, text = "Rule 2", variable = self._var_regle_2, \
                                onvalue = 1, offvalue = 0)
        C.grid(row=0, column=1)
        C = tk.Checkbutton(cd, text = "Rule 3", variable = self._var_regle_3, \
                                onvalue = 1, offvalue = 0)
        C.grid(row=0, column=2)
        C = tk.Checkbutton(cd, text = "Rule 4", variable = self._var_regle_4, \
                                onvalue = 1, offvalue = 0)
        C.grid(row=0, column=3)


        cd = tk.LabelFrame(cdr, text="Regressor:")
        cd.grid(row=0, column=1)
        self._var_regression = tk.IntVar(self._fn, value =0)

        R1 = tk.Radiobutton(cd, text="SVR", variable=self._var_regression, value=0)
        R1.grid(row=0, column=0)

        R2 = tk.Radiobutton(cd, text="KNN", variable=self._var_regression, value=1)
        R2.grid(row=0, column=1)

    def _calcul(self):
        self._L_mesures = list(self._list_mesures.get('@1,0', tk.END))
        if len(self._L_mesures)>0:
            self._L_ind = [0, len(self._L_mesures)]

            # Limites : 
            cond = True
            IC = float(self._var_IC.get())
            if IC > 99.999:
                IC = 99.999
            elif IC < 0.001:
                IC = 0.001
            M = int(self._var_nb_bootstrap.get())
            if M < 1:
                M = 1

            while cond:
                cond = False
                L_ind = self._L_ind
                for i in range(1, len(self._L_ind)):
                    ind_d = self._L_ind[i-1]
                    ind_f = self._L_ind[i]

                    id = bootstrap_cusum_p(self._L_mesures[ind_d:ind_f], k=0, M=M, 
                                      IC = IC,
                                      estimateur=int(self._var_esti_bootstrap.get()))
                    if id >= 0 and not (id+ind_d in L_ind):
                            cond = True
                            L_ind.append(id+ind_d)

                    id = bootstrap_cusum_m(self._L_mesures[ind_d:ind_f], k=0, M=M, 
                                      IC = IC,
                                      estimateur=int(self._var_esti_bootstrap.get()))
                    if id >= 0 and not (id+ind_d in L_ind):
                        cond = True
                        L_ind.append(id+ind_d)
                self._L_ind = sorted(L_ind)
            
            L_X_f = []
            L_Y_f = []
            L_X_m = []
            L_Y_m = []

            top = tk.Toplevel()
            plt.figure()
            for i in range(1,len(self._L_ind)):
                cdr = tk.LabelFrame(top, text="Group "+str(i)+ ":")
                print("Group "+str(i)+ ":")
                cdr.grid(row=i-1, column=0) 
                ind_d = self._L_ind[i-1]
                ind_f = self._L_ind[i]
                
                X = [i for i in range(ind_d, ind_f)]
                Y = self._L_mesures[ind_d:ind_f]

                Y_p = cusum_p(self._L_mesures[ind_d:ind_f])
                Y_m = cusum_m(self._L_mesures[ind_d:ind_f])
                plt.plot(X, Y_p[1:], '--b')
                plt.plot(X, Y_m[1:], '--g')
                lab = ttk.Label(cdr, text="Measurement index range: "+str(ind_d)+"-"+str(ind_f-1))
                print("Measurement index range: "+str(ind_d)+"-"+str(ind_f-1))
                lab.grid(row=0, column=0)
                if ind_f-ind_d>1:
                    mu = np.mean(self._L_mesures[ind_d:ind_f])
                    sigma = np.std(self._L_mesures[ind_d:ind_f])
                    print("Distribution : (µ, sigma)", mu, " ; ", sigma)
                    IC = int(10*float(self._var_monte_carlo_IC.get()))
                    if IC > 999:
                        IC = 999
                    elif IC < 1:
                        IC = 1
                    nb = int(self._var_nb_monte_carlo.get())
                    if nb<1:
                        nb = 1
                    taille = int(self._var_nb_taille_monte_carlo.get())
                    if taille < 1:
                        taille = 1
                    USL, LSL = monte_carlo(mu, sigma, 0, 
                                           IC=IC,
                                           nb=nb, 
                                           taille=taille)
                    lab = ttk.Label(cdr, text="Lower Control Limit: "+str(round(LSL,3)))
                    print("Lower Control Limit: "+str(round(LSL,3)))
                    lab.grid(row=1, column=0)
                    lab = ttk.Label(cdr, text="Upper Control Limit: "+str(round(USL,3)))
                    print("Upper Control Limit: "+str(round(USL,3)))
                    lab.grid(row=2, column=0)
                    plt.plot([min(X)-0.5, min(X)-0.5], [LSL, USL], '--k')
                    plt.plot([max(X)+0.5, max(X)+0.5], [LSL, USL], '--k')
                    plt.plot([min(X)-0.5, max(X)+0.5], [0, 0], '-k')
                    plt.plot([min(X)-0.5, max(X)+0.5], [USL, USL], '-r')
                    plt.plot([min(X)-0.5, max(X)+0.5], [USL*2/3, USL*2/3], '-', color=(1, 127/255, 39/255))
                    plt.plot([min(X)-0.5, max(X)+0.5], [USL/3, USL/3], '-', color=(1, 127/255, 39/255))
                    plt.plot([min(X)-0.5, max(X)+0.5], [LSL, LSL], '-r')
                    plt.plot([min(X)-0.5, max(X)+0.5], [LSL*2/3, LSL*2/3], '-', color=(1, 127/255, 39/255))
                    plt.plot([min(X)-0.5, max(X)+0.5], [LSL/3, LSL/3], '-', color=(1, 127/255, 39/255))
                    
                    X_err, Y_err = self._regles_p(X, Y_p[1:], USL)
                    L_ind = [ind for ind in range(len(X)) if not X[ind] in X_err]
                    L_X_reg = [X[ind] for ind in L_ind]
                    L_Y_reg = [Y[ind] for ind in L_ind]
                    L_X_err_reg = X_err
                    plt.plot(X_err, Y_err, 'or')
                    X_err, Y_err = self._regles_m(X, Y_m[1:], LSL)
                    L_ind = [ind for ind in range(len(X)) if not X[ind] in X_err and not X[ind] in L_X_reg]
                    L_X_reg += [X[ind] for ind in L_ind]
                    L_Y_reg += [Y[ind] for ind in L_ind]
                    L_X_err_reg += X_err
                    plt.plot(X_err, Y_err, 'or')
                    
                    L_X_m += X
                    L_Y_m += Y
                    if len(L_X_reg) > 2 and len(L_X_err_reg):
                        if int(self._var_regression.get())==0:
                            svr_rbf = SVR(kernel='rbf', C=100, gamma=0.1)
                            svr_rbf.fit(np.array(L_X_reg).reshape(-1,1), np.array(L_Y_reg).reshape(-1,1))
                            Y_reg = svr_rbf.predict(np.array(L_X_err_reg).reshape(-1,1))
                            
                            for j in range(len(L_X_err_reg)):
                                k = 0
                                while k < len(X):
                                    if X[k] == L_X_err_reg[j]:
                                        Y[k] = Y_reg[j]
                                        k = len(X)
                                    else:
                                        k += 1
                        else:
                            if len(L_X_reg)//3 > 1:
                                knn_regressor = KNeighborsRegressor(n_neighbors=len(L_X_reg)//3)
                            else:
                                knn_regressor = KNeighborsRegressor(n_neighbors=1)
                            knn_regressor.fit(np.array(L_X_reg).reshape(-1,1), np.array(L_Y_reg).reshape(-1,1))
                            Y_reg = knn_regressor.predict(np.array(L_X_err_reg).reshape(-1,1))
                            for j in range(len(L_X_err_reg)):
                                k = 0
                                while k < len(X):
                                    if X[k] == L_X_err_reg[j]:
                                        Y[k] = Y_reg[j][0]
                                        k = len(X)
                                    else:
                                        k += 1
                    L_X_f += X
                    L_Y_f += Y
            plt.xlabel("Time (s)")
            plt.ylabel("Cumulative Sum")
            plt.grid(axis="both")

            plt.figure()
            plt.plot(L_X_f, L_Y_f, '-xr', label='Change')
            plt.plot(L_X_m, L_Y_m, '-xk', label='Initial')
            plt.legend()
            plt.xlabel("Time (u.a.)")
            plt.ylabel("Value (u.a.)")

            plt.show()
            top.mainloop()

    def _regles_p(self, X, Y, L):
        X_f = []
        Y_f = []
        L_X3 = []
        L_Y3 = []
        L_X5 = []
        L_Y5 = []
        L_X8 = []
        L_Y8 = []
        for i in range(len(Y)):
            if Y[i] > L and int(self._var_regle_1.get())>0:
                X_f.append(X[i])
                Y_f.append(Y[i])
            if Y[i] > L*2/3 and int(self._var_regle_2.get())>0:
                L_X3.append(X[i])
                L_Y3.append(Y[i])
                if len(L_X3)>2:
                    L_X3.pop(0)
                    L_Y3.pop(0)
                if len(L_X3) >=2:
                    X_f += L_X3
                    Y_f += L_Y3
            else:
                if len(L_X3)>0:
                    L_X3.pop(0)
                    L_Y3.pop(0)
            if Y[i] > L/3 and int(self._var_regle_3.get())>0:
                L_X5.append(X[i])
                L_Y5.append(Y[i])
                if len(L_X5)>4:
                    L_X5.pop(0)
                    L_Y5.pop(0)
                if len(L_X5) >=4:
                    X_f += L_X5
                    Y_f += L_Y5
            else:
                if len(L_X5)>0:
                    L_X5.pop(0)
                    L_Y5.pop(0)

            if Y[i] > L/10 and int(self._var_regle_4.get())>0:
                L_X8.append(X[i])
                L_Y8.append(Y[i])
                if len(L_X8)>7:
                    L_X8.pop(0)
                    L_Y8.pop(0)
                if len(L_X8) >=7:
                    X_f += L_X8
                    Y_f += L_Y8
            else:
                if len(L_X8)>0:
                    L_X8.pop(0)
                    L_Y8.pop(0)
        return X_f, Y_f

    def _regles_m(self, X, Y, L):
        X_f = []
        Y_f = []
        L_X3 = []
        L_Y3 = []
        L_X5 = []
        L_Y5 = []
        L_X8 = []
        L_Y8 = []
        for i in range(len(Y)):
            if Y[i] < L and int(self._var_regle_1.get())>0:
                X_f.append(X[i])
                Y_f.append(Y[i])
            if Y[i] < L*2/3 and int(self._var_regle_2.get())>0:
                L_X3.append(X[i])
                L_Y3.append(Y[i])
                if len(L_X3)>2:
                    L_X3.pop(0)
                    L_Y3.pop(0)
                if len(L_X3) >=2:
                    for elmt in L_X3:
                        if not elmt in X_f:
                            X_f.append(elmt)
                            Y_f.append(elmt) 
            else:
                if len(L_X3)>0:
                    L_X3.pop(0)
                    L_Y3.pop(0)
            if Y[i] < L/3 and int(self._var_regle_3.get())>0:
                L_X5.append(X[i])
                L_Y5.append(Y[i])
                if len(L_X5)>4:
                    L_X5.pop(0)
                    L_Y5.pop(0)
                if len(L_X5) >=4:
                    X_f += L_X5
                    Y_f += L_Y5
            else:
                if len(L_X5)>0:
                    L_X5.pop(0)
                    L_Y5.pop(0)

            if Y[i] < L/10 and int(self._var_regle_4.get())>0:
                L_X8.append(X[i])
                L_Y8.append(Y[i])
                if len(L_X8)>7:
                    L_X8.pop(0)
                    L_Y8.pop(0)
                if len(L_X8) >=7:
                    X_f += L_X8
                    Y_f += L_Y8
            else:
                if len(L_X8)>0:
                    L_X8.pop(0)
                    L_Y8.pop(0)
        return X_f, Y_f
    
if '__main__' ==__name__:
    fn = app()  