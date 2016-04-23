### build.jl --- Build script for CmplxRoots.jl.

# Copyright (C) 2016  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

### Code:

local_dir  = dirname(@__FILE__)
local_file = "cmplx_roots_sg_77.f"
libcmplxroots = joinpath(local_dir, "libcmplxroots")
object=""; opts=""
@linux_only (object="libcmplxroots.so";    opts="-shared")
@osx_only   (object="libcmplxroots.dylib"; opts="-dynamiclib")

# Download cmplx_roots_sg and build the shared object.
info("Downloading cmplx_roots_sg source...")
download("http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/cmplx_roots_sg_77.f",
         local_file)
info("Building libcmplxroots...")
cd(local_dir) do
    run(`gfortran $opts -fPIC -o $object $local_file`)
end

# Make sure Julia is able to see the library.
if length(Libdl.find_library([libcmplxroots])) > 0
    info("libcmplxroots successfully installed!")
else
    error("Installation of libcmplxroots failed")
end
