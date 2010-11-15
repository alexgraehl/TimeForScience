all: $(DATA_FILE) $(GXA_FILE)

clean:
	$(RM) $(GXA_FILE) $(DATA_FILE) $(SQL_DSC_FILE) $(SQL_DATA_FILE)

remote:
	$(TYPICAL_REMOTE)

remote_clean:
	$(RM) $(REMOTE_DIR)/*

$(DATA_FILE): $(REMOTE_TARGET)
	cd ../; parse_go.pl $(ORGANISM) $(ORGANISM)/$(REMOTE_TARGET); \
	cat $(ORGANISM).go.dat | \
	$(UPPERCASE_KEYS) | \
	$(COMMON_GO_PIPES) > \
	$(ORGANISM)/$(DATA_FILE); \
	$(RM) $(FLAT_GO_FILE); \
	cd $(ORGANISM); \

$(GXA_FILE): $(DATA_FILE)
	tab2gxa.pl $(DATA_FILE) > $(GXA_FILE)

