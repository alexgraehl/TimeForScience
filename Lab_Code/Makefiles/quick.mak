default: all


all:
	@echo "You should use recurse.mak, and not quick.mak. It is the same in almost all respects, but recurse.mak has much better error handling."
	$(error "USE RECURSE.MAK, NOT QUICK.MAK!!!!")
