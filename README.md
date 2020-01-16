# EZ_PARALLEL: A Fortran 95 Library for Simplifying Parallelizing 2D Finite Difference Codes

## Background

EZ_PARALLEL is a library for Fortran 95 created to help simplify the processes of parallelizing serial finite difference codes, by including pre-made subroutines for grid decomposition, adjacent subgrid communication, and more. We cannot guarantee that this library will yield the most optimal results, but our hope is that people without the time to learn and/or experience with MPI can use this library to improve the scope of their work.

EZ_PARALLEL was developed by Jason Turner under the guidance of Samuel Stechmann at the University of Wisconsin-Madison. 

## How Does EZ_PARALLEL Work?

EZ\_PARALLEL acts as an interface between the user and MPI to make the user's experience with MPI simpler. As with any interface of this nature, it cannot do everything a user might hope for. When writing a parallel code with MPI, each processor will execute the entire code. This can present an issue in some areas of serial codes, such as outputting data since the output file name will need to be different for each processor. (Otherwise, each processor will keep overwriting the same output files!) Although we aim to make this library as general as possible, there are some ways in which the serial code must function in order for this library to be used. 

Generally, we expect the flow of a serial finite difference code to be as follows:
1. Initialize simulations parameters, e.g., the size of the array containing the domain (with allocating memory to that array), setting physical constants, etc.
1. Initializing the domain grid. This includes allocating memory to the array containing the domain and setting the initial condition.
1. Step forward in time. This includes calculating each subsequent time step and outputting the domain grid to a file.

This is not to say that every code must be in the format to use EZ\_PARALLEL, but it is the basic form that we expect any basic serial code to be in. We will explain the more specific requirements in the context of the relevant EZ\_PARALLEL subroutines below:

### Memory Allocation for the Domain Array

The GRID\_DECOMP subroutine takes in the domain size, say (x\_len, y\_len), and assigns to those values a subdomain size based on the original domain size, number of processors, and the number of grid points each subdomain would need to borrow from its adjacent subdomains. Note, each subdomain retains the same x\_len value as the domain, which may present difficulties for very skewed domains.

We shall explain how GRID\_DECOMPOSITION works using a simple example. The requirement for the serial code is stated explicitly at the end of this subsection. Let's say we have a domain of size (total\_x\_len, total\_y\_len) = (9, 11) to be divided among three processors: processor 0, processor 1, and processor 2. First, we give each processor a y\_len of FLOOR(total\_y\_len/proc\_count), where proc\_count is the total number of processors. However, these subdomains do not cover the entire domain (each processor has y\_len = 3).

To resolve this, we add 1 to the y\_len of each processor with ID less than MODULO(total\_y\_len, proc\_count). In this case, processors 0 and 1 have y\_len = 4, while processor 2 still has y\_len = 3. Geometrically, we consider processor 0 to contain the bottom-most rows and the final processor (the one with the greatest ID, in this case 2) to contain the top-most rows. Processor 1 will be immediately above processor 0, processor 2 will be immediately above processor 1, and so forth.

Now, let U<sub>i,j</sub><sup>n+1</sup> denote the numerical solution at time t<sub>n+1</sub> and point (x<sub>i</sub>, y<sub>j</sub>). Any numerical scheme will require points U<sub>i+p,j+q</sub><sup>n+k</sup> to calculate U<sub>i,j</sub><sup>n+1</sup> (p, q, k being integers). The subgrids may not contain the necessary points to calculate the value of U<sub>i,j</sub><sup>n+1</sup> at each point in the subgrid, and may need to borrow some points from its neighbors. For our example, let's say we may calculate U<sub>i,j</sub><sup>n+1</sup> using U<sub>i+1,j</sub><sup>n</sup>, U<sub>i-1,j</sub><sup>n</sup>, U<sub>i,j+1</sub><sup>n</sup>, U<sub>i,j-1</sub><sup>n</sup>, and U<sub>i,j</sub><sup>n</sup>.

Processor 0 contains the points (x<sub>i</sub>, y<sub>j</sub>) for 1 &le; i &le; 9, 1 &le; j &le; 4. To calculate U<sub>i,4</sub><sup>n</sup>, processor 0 will need the value of the numerical solution at (x<sub>i</sub>, y<sub>5</sub>), which is stored on processor 1. To address this, we add 1 to the value of y\_len of processor 0 to store the values U<sub>i,5</sub><sup>n</sup> obtained from processor 1. Likewise, we add 1 to the value of y\_len of processor 2 to store the values of U<sub>i,8</sub><sup>n</sup> from processor 1.

We will need to add 2 to the value of y\_len for processor 1, since it must store the values U<sub>i,4</sub><sup>n</sup> from processor 0 and U<sub>i,9</sub><sup>n</sup> from processor 2.

More generally, we say that this numerical scheme has an _overlap_ of 1. For exterior subgrids (that with processor ID 0 and that with the largest processor ID), we must add the overlap value to y\_len to accomodate for the points they must borrow from adjacent subdomains. For interior subgrids (all others), we must add twice the overlap value to y\_len to accomodate for the points they must borrow from their neighbors. We call the points each processor must borrow the _overlap points for the subdomain of that processor_, and all other points in the subdomain of the processor the _interior points for the subdomain of that processor_.

To complete the example, processor 0 has y\_len = 5, processor 1 has y\_len = 6, and processor 2 has y\_len = 4.

Therefore, the array used to store the domain grid in the serial code must be allocatable, and must be allocated using two parameters. The two parameters will be adjusted by GRID\_DECOMP so that the appropriate-sized arrays are allocated for each processor.

### Setting the Initial Condition and Dirichlet Boundary Conditions

There are many ways to set an initial condition in a code: using the values from a specific function, reading a state from an input file, generating random noise, and more. Since each processor in a parallel code attempts to execute the entire code, the serial subroutines for setting the initial condition generally do not translate well to each subdomain.

Unfortunately, EZ\_PARALLEL does not contain a subroutine for reading an initial condition from a file. However, it can support setting the initial condition for each subdomain using an initial condition function.

The serial code must have a reference point (located in the bottom-left corner of the domain) to set its initial condition, as well as uniform grid spacing. The INDENTIFY\_REF\_POINT subroutine of the EZ\_PARALLEL library will identify the location of the reference point in space for each subdomain, and change that processor's value of (x\_ref, y\_ref) accordingly.

Continuiung from our previous example, say this refence point is (x\_ref, y\_ref) = (0, 0), the grid spacing is dx in the x-direction and dy in the y-direction. Since the bottom-left point of the subdomain for processor 0 is also the bottom-left point of the domain, processor 0 has (x\_ref, y\_ref). The bottom-left point of the subdomain for processor 1 corresponds to the point (x<sub>0</sub>, y<sub>4</sub>) = (0, 4\*dy). Thus, we set (x\_ref, y\_ref) = (0, 4\*dy) for processor 1. To complete the example, we set (x\_ref, y\_ref) = (0, 8\*dy) for processor 2.

Setting Dirichlet boundary conditions provides a very similiar difficulty to setting the initial condition. Since each processor calls the entire code, the boundary of each subdomain will have Dirichlet boundary conditions (meaning that interior points of the domain are assigned Dirichlet conditions!). However, since the boundary of each subdomain contained in the overlap points, we may simply call SHARE\_SUBGRID\_BOUNDARIES after each time we set the Dirichlet boundary conditions to ammend this issue.