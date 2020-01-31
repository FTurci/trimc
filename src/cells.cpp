
#include "cells.h"

CellList::CellList() : dimension(3)
{
}

CellList::CellList(unsigned int dimension_, const std::vector<double>& boxSize, const std::vector<double>& boxCentre, double range) : dimension(dimension_)
{
    this->initialise(boxSize, centre,range);
}

CellList& CellList::operator = (const CellList& cells)
{
    // resize cell list
    (*this).resize(cells.size());

    // copy across individual cells
    for (unsigned int i=0;i<cells.size();i++)
    {
        at(i) = cells[i];
    }

    return *this;
}

CellList& CellList::operator = (const std::vector<Cell>& cells)
{
    // resize cell list
    (*this).resize(cells.size());

    // copy across individual cells
    for (unsigned int i=0;i<cells.size();i++)
    {
        at(i) = cells[i];
    }

    return *this;
}

void CellList::initialise(const std::vector<double>& boxSize, const std::vector<double>& boxCentre, double range)
{
    unsigned int i,j,k,m;
    unsigned int a,b,c;
    unsigned int x,y,z;
    unsigned int nn,nnCount;

    cellsPerAxis.resize(dimension);
    cellSpacing.resize(dimension);
    centre.resize(dimension);

    // printf("dim %d size %g range %g \n",dimension, boxSize[0], range);

    for (i=0;i<dimension;i++)
    {   
        // centre of the cell system: if the box has centre in {0,0,0} the cell system has centre in {L/2 L/2 L/2}
        centre[i] = +boxSize[i]/2-boxCentre[i] ;
        // printf ("#dimension %d size %g\n", i ,boxSize[i]);
        cellsPerAxis[i] = 1;

        while ((boxSize[i] / (double) cellsPerAxis[i]) > range)
        {
            cellsPerAxis[i]++;
        }
        cellsPerAxis[i]--;
        cellSpacing[i] = boxSize[i] / (double) cellsPerAxis[i];

        // check that number of cells per axis is large enough
        if (cellsPerAxis[i] < 3)
        {
            std::cerr << "[ERROR] CellList: Simulation box is too small (min cells per axis is 3)\n";
            exit(EXIT_FAILURE);
        }
        printf(":: Axis %d : %d cells\n", i,cellsPerAxis[i]);
    }

    // Estimate maximum number of particles per cell from interaction range.
    // (Assumes particle diameter is one.)
    if (dimension == 3) maxParticles = (cellSpacing[0]*cellSpacing[1]*cellSpacing[2]) / ((4.0/3.0)*M_PI*0.5*0.5*0.5);
    else maxParticles = (cellSpacing[0]*cellSpacing[1]) / (M_PI*0.5*0.5);

    // Add a buffer, e.g. if particles can overlap.
    maxParticles += 100;
    printf(":: Max number of particles %d.\n", maxParticles);
    nCells = cellsPerAxis[0]*cellsPerAxis[1];
    if (dimension == 3) nCells *= cellsPerAxis[2];

    // resize cell list array
    (*this).resize(nCells);

    if (dimension == 3)
    {
        nNeighbours = 27;

        // loop over all cells x direction
        for (i=0;i<cellsPerAxis[0];i++)
        {
            // loop over all cells y direction
            for (j=0;j<cellsPerAxis[1];j++)
            {
                // loop over all cells z direction
                for (k=0;k<cellsPerAxis[2];k++)
                {
                    // cell index
                    m = i + cellsPerAxis[0]*j + cellsPerAxis[0]*cellsPerAxis[1]*k;

                    // resize neighbour list and particle arrays
                    at(m).neighbours.resize(nNeighbours);
                    at(m).particles.resize(maxParticles);

                    nnCount = 0;

                    // x loop for nearest neighbours
                    for (a=0;a<3;a++)
                    {
                        x = (i+(a-1)+cellsPerAxis[0])%cellsPerAxis[0];

                        // y loop for nearest neighbours
                        for (b=0;b<3;b++)
                        {
                            y = (j+(b-1)+cellsPerAxis[1])%cellsPerAxis[1];

                            // z loop for nearest neighbours
                            for (c=0;c<3;c++)
                            {
                                z = (k+(c-1)+cellsPerAxis[2])%cellsPerAxis[2];

                                // nn cell index
                                nn = x + y*cellsPerAxis[0] + z*cellsPerAxis[0]*cellsPerAxis[1];

                                at(m).neighbours[nnCount] = nn;
                                nnCount++;
                            }
                        }
                    }

                    // zero tally for each cell
                    at(m).tally = 0;
                    at(m).index = m;
                }
            }
        }
    }
    else
    {
        nNeighbours = 9;

        // loop over all cells x direction
        for (i=0;i<cellsPerAxis[0];i++)
        {
            // loop over all cells y direction
            for (j=0;j<cellsPerAxis[1];j++)
            {
                // cell index
                m = i + cellsPerAxis[0]*j;

                // resize neighbour list and particle arrays
                at(m).neighbours.resize(nNeighbours);
                at(m).particles.resize(maxParticles);

                nnCount = 0;

                // x loop for nearest neighbours
                for (a=0;a<3;a++)
                {
                    x = (i+(a-1)+cellsPerAxis[0])%cellsPerAxis[0];

                    // y loop for nearest neighbours
                    for (b=0;b<3;b++)
                    {
                        y = (j+(b-1)+cellsPerAxis[1])%cellsPerAxis[1];

                        // nn cell index
                        nn = x + y*cellsPerAxis[0];

                        at(m).neighbours[nnCount] = nn;
                        nnCount++;
                    }
                }

                // zero tally for each cell
                at(m).tally = 0;
                at(m).index = m;
            }
        }
    }
}

void CellList::reset()
{
    for (unsigned int i=0;i<nCells;i++) at(i).tally = 0;
}

int CellList::getCell(const Particle& particle)
{
    int cell,cellx,celly;

    double x = particle.pos[0]+centre[0];
    double y = particle.pos[1]+centre[1];
    double z ;
    if (dimension == 3)
        z = particle.pos[2]+centre[2];
    // printf("x y z %g %g %g\n",x,y,z);

    cellx = int(x/cellSpacing[0]);
    celly = int(y/cellSpacing[1]);

    cell = cellx + celly*cellsPerAxis[0];
    // printf("cell %d %d %d\n",cellx,celly,cell);
    if (dimension == 3)
    {
        int cellz = int(z/cellSpacing[2]);
        cell += cellz*cellsPerAxis[0]*cellsPerAxis[1];


    }

    return cell;
}

void CellList::initCell(int newCell, Particle& particle)
{
    // Add to new list
    // printf ("This particle index %d newcell %d\n", particle.index,newCell);
    at(newCell).particles[at(newCell).tally] = particle.index;
    // printf("Accessed\n");
    particle.cell = newCell;
    particle.posCell = at(newCell).tally;
    at(newCell).tally++;

    if (at(newCell).tally == maxParticles)
    {
        std::cerr << "[ERROR] CellList: Maximum number of particles per cell exceeded!\n";
        exit(EXIT_FAILURE);
    }
}

void CellList::initCellList(std::vector<Particle>& particles)
{
    for (unsigned int i=0;i<particles.size();i++)
    {
        initCell(getCell(particles[i]), particles[i]);
    }
}

void CellList::updateCell(int newCell, Particle& particle, std::vector<Particle>& particles)
{
    // Remove from old list
    at(particle.cell).tally--;
    at(particle.cell).particles[particle.posCell] = at(particle.cell).particles[at(particle.cell).tally];
    particles[at(particle.cell).particles[at(particle.cell).tally]].posCell = particle.posCell;

    // Add to new list
    at(newCell).particles[at(newCell).tally] = particle.index;
    particle.cell = newCell;
    particle.posCell = at(newCell).tally;
    at(newCell).tally++;

    if (at(newCell).tally == maxParticles)
    {
        std::cerr << "[ERROR] CellList: Maximum number of particles per cell exceeded!\n";
        exit(EXIT_FAILURE);
    }
}

void CellList::setDimension(unsigned int dimension_)
{
    dimension = dimension_;
}

unsigned int CellList::getNeighbours() const
{
    return nNeighbours;
}