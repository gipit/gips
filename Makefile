install_dev:
	cp GIP/gip/bin/Release/libgip.so /usr/lib/	
	cp GIP/giputils/bin/Release/* /usr/local/bin/
	#cp GIP/gippy/bin/Release/_gippylib.so /usr/lib/python2.7/
	#cp GIP/gippy/gippylib.py /usr/lib/python2.7/
	#ln -s `pwd`/bin/* /usr/local/bin/
	#ln -s `pwd`/gippy /usr/lib/python2.7

install:
	cp GIP/gip/bin/Release/libgip.so /usr/lib/	
	cp GIP/giputils/bin/Release/* /usr/local/bin/
	cp GIP/gippy/bin/Release/_gippylib.so /usr/lib/python2.7/
	cp GIP/gippy/gippylib.py /usr/lib/python2.7/
	cp bin/* /usr/local/bin/
	cp -R gippy /usr/lib/python2.7/
	pip install -r requirements.txt

clean:
	rm /usr/lib/libgip.so
	rm /usr/local/bin/gip*
	rm /usr/lib/python2.7/_gippylib.so
	rm /usr/lib/python2.7/gippylib.py
	rm -rf /usr/lib/python2.7/gippy
	rm /usr/local/bin/landsat
