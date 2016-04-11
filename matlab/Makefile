all: genVelPerturb

genVelPerturb : genVelPerturb.m
	mcc -m $<

clean:
	rm -f readme.txt mcc*.log  run_*.sh genVelPerturb
