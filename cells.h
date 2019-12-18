// Copied and adapted from Lester Hedges's VMMC lester.hedges+vmmc@gmail.com>

#ifndef _CELLS_H
#define _CELLS_H


#include <iostream>
#include <vector>
#include "particle.h"


//! Structure containing attributes for an individual cell.
struct Cell
{
    unsigned int index;                     //!< Cell index.
    unsigned int tally;                     //!< Number of particles in the cell.
    std::vector<unsigned int> particles;    //!< Indices of particles in the cell.
    std::vector<unsigned int> neighbours;   //!< Indices of nearest neighbour cells.
};

//! Container class for storing a list of cells.
//! This class contains the main cell list that is manipulated by the simulation.
class CellList : public std::vector<Cell>
{
public:
    //! Default constructor.
    CellList();

    //! Constructor.
    /*! \param dimension_
            The dimensionality of the simulation box.

        \param boxSize
            The size of the simulation box in each dimension.

        \param range
            Maximum interaction range.
     */
    CellList(unsigned int, const std::vector<double>&,const std::vector<double>& , double);

    //! Copy constructor. \param cells A reference to an existing CellList object.
    CellList& operator = (const CellList&);

    //! Copy constructor. \param cells A vector of existing Cell data structures.
    CellList& operator = (const std::vector<Cell>&);

    //! Initialise cell lists.
    /*! \param boxSize
            The size of the simulation box in each dimension.

        \param range
            Maximum interaction range.
     */
    void initialise(const std::vector<double>&, const std::vector<double>& ,double );

    //! Reset cell lists (zero cell tallys).
    void reset();

    //! Get cell index for a particle.
    /*! \param particle
            Reference to a particle.

        \return
            The cell index.
     */
    int getCell(const Particle&);

    //! Initialise cell list for an individual particle.
    /*! \param newCell
            The index of the cell in which the particle is located.

        \param particle
            Reference to a particle.
     */
    void initCell(int, Particle&);

    //! Initialise cell list for all particles.
    /*! \param particles Reference to a vector of particles.
     */
    void initCellList(std::vector<Particle>&);

    //! Update cell list for an individual particle.
    /*! \param newCell
            The index of the cell in which the particle is located.

        \param particle
            Reference to a particle.

        \param particles
            Reference to a vector of particles.
     */
    void updateCell(int, Particle&, std::vector<Particle>&);

    //! Set the dimensionality of the cell list.
    /*! \param dimension_
            The dimensionality of the simulation.
     */
    void setDimension(unsigned int);

    //! Get the number of neighbours per cell.
    unsigned int getNeighbours() const;

private:
    std::vector<double> centre;                 //!< Centre of the simulation box
    unsigned int dimension;                     //!< Dimension of the simulation box.
    unsigned int nCells;                        //!< Total number of cells.
    unsigned int nNeighbours;                   //!< Number of neighbours per cell.
    unsigned int maxParticles;                  //!< Maximum number of particles per cell.
    std::vector<unsigned int> cellsPerAxis;     //!< Number of cells per axis.
    std::vector<double> cellSpacing;            //!< Spacing between cells.
};

#endif