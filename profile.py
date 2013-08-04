import cProfile
command = """import main"""
cProfile.runctx( command, globals(),locals(),filename="aestimo_numpy-0.8.3.profile")

#command = """import aestimo"""
#cProfile.runctx( command, globals(), locals(), filename="aestimo_t8.profile" )

#command = """import main"""
#cProfile.runctx( command, globals(), locals(), filename="aestimo_numpy.profile" )