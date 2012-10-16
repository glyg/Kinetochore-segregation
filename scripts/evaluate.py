from kt_simul.analysis import evaluations as eva

evals = eva.find_evaluations()
print evals

rate = evals[0]()
rate.run()