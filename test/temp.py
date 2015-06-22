from navierstokes_paper1 import *
state = StateEqn(2,2,np.array([0,1]),np.array([0,1]))
coef = array([[0,1],[.5,.3]])
state.setcoef(coef)
U = array([1,2,1])
R = array([3,4,5])
Pres = state.P(U,R)[0]

Ures = state.U(Pres,R)
