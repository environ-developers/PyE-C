test: test-serial

test-serial:
	@ mpirun -n 1 python3 tester.py

test-parallel:
	@ mpirun -n 4 python3 tester.py

clean:
	@ rm environ.debug CRASH *.out