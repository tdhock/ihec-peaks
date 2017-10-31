test.out: test.bedGraph denormalize
	./denormalize test.bedGraph test.out
denormalize: denormalize.cpp
	g++ denormalize.cpp -o denormalize
