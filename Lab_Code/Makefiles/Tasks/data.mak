
ifndef ($(targets))
   targets = data.tab
endif

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)

# To be effectual, need to define data.tab (i.e. fill
# in the data.tab rule at the bottom).

ifneq (,$(REMOTE_INFO))

remote:
	$(MKDIR) $(REMOTE_DIR); \
	$(foreach u, $(URL), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) '$(u)' $(URL_INFO); \
	) \
	$(RM) $(REMOTE_DIR)/$(README); \
	cd $(REMOTE_DIR); \
	$(LINK) ../$(REMOTE_INFO) $(README); \

else

remote:
	$(MKDIR) $(REMOTE_DIR); \
	$(foreach u, $(URL), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) '$(u)'; \
	) \

endif

remote_clean:
	$(RM) $(REMOTE_DIR)/*

make:

