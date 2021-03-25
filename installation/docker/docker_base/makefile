# A MAKE FILE TO COPY THE RIGHT CODES AND CREATE THE DOCKER IMAGE
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)

name = "fabiankindermann/ce-fortran"

# CLEAN UP THE WORKSPACE
.PHONY: clean
clean :
	rm -rf code-book
	rm -rf code-solution
	rm -rf toolbox.f90

.PHONY: build
build :
	cp -R ../../../code-book .
	cp -R ../../../code-solution .
	cp ../../toolbox/toolbox.f90 .
	docker build -t $(name) . 
	docker tag $(name) $(name):$(version)
	rm -rf code-book
	rm -rf code-solution
	rm -rf toolbox.f90