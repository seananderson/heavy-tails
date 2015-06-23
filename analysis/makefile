CACHE = ~/Dropbox/heavy-tails/cache

all: analysis

analysis:
	Rscript --vanilla make.R

# Warning: this will take a *long* time to rebuild after:
clean:
	rm *.rds
	rm *.rda

sendCacheToDropbox:
	cp gpdd-clean.rds $(CACHE)
	cp gomp-base-hat.rds $(CACHE)
	cp logistic-hat.rds $(CACHE)
	cp gomp-ar1-hat.rds $(CACHE)
	cp gomp-obs-0.2-hat.rds $(CACHE)
	cp rate-hat.rds $(CACHE)
	cp gomp-base-weaker-hat.rds $(CACHE)
	cp gomp-base-stronger-hat.rds $(CACHE)
	cp gomp-base-mean-sd.rds $(CACHE)
	cp beta-stan-samples.rds $(CACHE)
	cp beta-stan-samples-n-only.rds $(CACHE)
