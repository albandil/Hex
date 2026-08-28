all: bin/ILUP

bin/ILUP: src/ILUP.cpp
	@mkdir -p $(@D)
	g++ $< -o $@

clean:
	rm -f bin/ILUP
