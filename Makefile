install:
	cp GIP/gip/bin/Release/libgip.so /usr/lib/
	#cp GIP/giputils/bin/Release/* /usr/local/bin/
	#git clean -xfd
	#python setup.py install

develop:
	cp GIP/gip/bin/Release/libgip.so /usr/lib/
	git clean -xfd
	python setup.py develop

clean:
	rm /usr/lib/libgip.so
	#rm /usr/local/bin/gip_*
