==============================================
GridPACK\ :math:`\mathrm{{}^{TM}}` User Manual
==============================================

:Author: B. Palmer, W. Perkins, S. Abhyankar
:Date: 2025-01-10

**Abstract**: This document is designed to provide an in-depth
description of the GridPACK\ :math:`\mathrm{{}^{TM}}` framework and the
software modules contained within it. In combination with the
Doxygen-based documentation on the GridPACK webpage, users and
application developers should have a complete description of the
framework components and how to use them. The applications area in the
source code directory, as well as the GridPACK modules and components,
provide additional examples of how GridPACK can be used to create power
grid applications. However, if there are still questions on GridPACK,
users should feel free to contact the GridPACK development team.

**License**

GridPACK\ :math:`\mathrm{{}^{TM}}` License: Copyright (c) 2013, Battelle
Memorial Institute All rights reserved.

1. Battelle Memorial Institute (hereinafter Battelle) hereby grants
permission to any person or entity lawfully obtaining a copy of this
software and associated documentation files (hereinafter "the Software")
to redistribute and use the Software in source and binary forms, with or
without modification. Such person or entity may use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the
Software, and may permit others to do so, subject to the following
conditions: \* Redistributions of source code must retain the above
copyright notice, this list of conditions and the following disclaimers.
\* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution. \*
Other than as used herein, neither the name Battelle Memorial Institute
or Battelle may be used in any form whatsoever without the express
written consent of Battelle.

2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

3. The Software was produced by Battelle under Contract No.
DE-AC05-76RL01830 with the Department of Energy. For five years from
October 10, 2013, the Government is granted for itself and others acting
on its behalf a nonexclusive, paid-up, irrevocable worldwide license in
this data to reproduce, prepare derivative works, and perform publicly
and display publicly, by or on behalf of the Government. There is
provision for the possible extension of the term of this license.
Subsequent to that period or any extension granted, the Government is
granted for itself and others acting on its behalf a nonexclusive,
paid-up, irrevocable worldwide license in this data to reproduce,
prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so. The
specific term of the license can be identified by inquiry made to
Battelle or DOE. Neither the United States nor the United States
Department of Energy, nor any of their employees, makes any warranty,
express or implied, or assumes any legal liability or responsibility for
the accuracy, completeness or usefulness of any data, apparatus, product
or process disclosed, or represents that its use would not infringe
privately owned rights.

**How to read this document**

Depending on how you are planning on using
GridPACK\ :math:`\mathrm{{}^{TM}}`, there are a variety of different
ways of approaching the documentation. If you are only planning on using
existing applications as is, without modification, then you should focus
on the sections for configuring and building GridPACK and the
application module documentation. Users that are interested in
developing their own applications may want to scan the section
“Developing Applications” before going to the beginning of the document
to learn about individual functionality in depth. The “GridPACK
Examples” section contains additional examples of simple applications
that can be used to get a sense of how to build an application from the
ground up. The “Contingency Analysis” section provides some information
on how to build applications that are based on the existing GridPACK
application modules. 
