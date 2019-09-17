#functions written to by me
import functions as fun
import Vegas_final as vg
import simpsons_trapezium_final as st
import transformation_final as tr

# to plot iterations below
import matplotlib.pyplot as plt

st.solve_ts() # solve simpsons trapezium for wavefunction (produces 2 plots)
input('Press any Key to continue (1)')

st.val_ts(plot = False) # validate simpsons trapezium for three functions
# st.val_ts(plot = True) to reproduce 6 plots
input('Press any Key to continue (2)')

# solve monte carlo importance sampling for uniform pdf to rel accuracy 1e-6 (produces 2 plots)
tr.solve_tran(pdf = fun.uni_pdf, trans=fun.trans_uni, e = 1e-6)
input('Press any Key to continue (3)')

# solve monte carlo importance sampling for linear pdf to rel accuracy 1e-6 (produces 2 plots)
tr.solve_tran(e=1e-6)
input('Press any Key to continue (4)')

# validate monte carlo importance sampling via tranfomration method
tr.val_tran(plot = False)
#tr.val_tran(plot = True) to reproduce 6 plots
input('Press any Key to continue (5)')

tr.plot_times() # 'plot the iterations taken with relative accuracy

#%%

# solve Vegas for relative accuracy of 1e-4 (produces 2 plots)
vg.solve_veg(e = 1e-4)
input('Press any Key to continue (7)')

# validate Vegas on lower than previous relative accuracies for reasonable
# run times
vg.val_veg(plot = False)
#vg.val_veg(plot = True) to reproduce 6 plots

