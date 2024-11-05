import IPython
program = "obpm.py"
process = 2
command = "!mpiexec -n %d python %s" % (process, program)
IPython.get_ipython().run_cell(command)
