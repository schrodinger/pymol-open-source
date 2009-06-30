/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: unit_conversion.h,v $
 *      $Author: akohlmey $       $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $       $Date: 2009/06/26 01:13:39 $
 *
 ***************************************************************************/
/******************************************************************
 * 
 * unit conversion factors.
 *  
 ******************************************************************/


#ifndef UNIT_CONVERSION_H
#define UNIT_CONVERSION_H

/* 
 * units according to CODATA 2006 
 * Version 5.1, http://physics.nist.gov/constants
 */
/* convert Bohr to Angstrom */
#define BOHR_TO_ANGS 0.529177210859
#define ANGS_TO_BOHR 1.88972612478289694072


/* convert Hartree into kcal/mol */
#define HARTREE_TO_KCAL 627.5094706142

/* convert Hartree into eV */
#define HARTREE_TO_EV    27.211383860484776

#endif /* UNIT_CONVERSION */

