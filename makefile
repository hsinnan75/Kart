.KEEP_STAT:

all: kart bwt_index

kart: htslib
		$(MAKE) -C src && cp -f src/$@ .

htslib:
		$(MAKE) -C src/htslib libhts.a

bwt_index:
		$(MAKE) -C src/BWT_Index && cp -f src/BWT_Index/$@ .

clean:
		rm -f kart bwt_index
		$(MAKE) clean -C src
		$(MAKE) clean -C src/htslib
		$(MAKE) clean -C src/BWT_Index
