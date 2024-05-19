import pinocchio as pin
import numpy as np
import time

model = pin.buildSampleModelManipulator()
data = model.createData()
 
q = pin.neutral(model)
v = pin.utils.zero(model.nv)
a = pin.utils.zero(model.nv)
 
tau = pin.rnea(model,data,q,v,a)
print('tau = ', tau.T)