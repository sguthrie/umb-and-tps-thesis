# -*- encoding: utf-8 -*-
# Hydrogen Bonding Analysis
"""
Hydrogen Bond analysis --- :mod:`MDAnalysis.analysis.hbonds`
===================================================================

:Author: David Caplan, Lukas Grossar, Oliver Beckstein
:Year: 2010-2012
:Copyright: GNU Public License v3


Given a :class:`~MDAnalysis.core.AtomGroup.Universe` (simulation
trajectory with 1 or more frames) measure all hydrogen bonds for each
frame between selections 1 and 2.

The :class:`HydrogenBondAnalysis` class is modeled after the `VMD
HBONDS plugin`_.

.. _`VMD HBONDS plugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/hbonds/

Options:
  - *update_selections* (``True``): update selections at each frame?
  - *selection_1_type* ("both"): selection 1 is the: "donor", "acceptor", "both"
  - donor-acceptor *distance* (Å): 3.0
  - Angle *cutoff* (degrees): 120.0
  - *forcefield* to switch between default values for different force fields
  - *donors* and *acceptors* atom types (to add additional atom names)

.. _Analysis Output:

Output
------

The results are hydrogen bond data per frame (# indicates comments
that are not part of the output), stored in
:attr:`HydrogenBondAnalysis.timeseries`::

    results = [
        [ # frame 1
           [ # hbond 1
              <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
           ],
           [ # hbond 2
              <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
           ],
           ....
        ],
        [ # frame 2
          [ ... ], [ ... ], ...
        ],
        ...
    ]

.. Note::

   For historic reasons, the *donor index* and *acceptor index* are a 1-based
   indices. To get the :attr:`Atom.number` (the 0-based index typically used in
   MDAnalysis simply subtract 1. For instance, to find an atom in
   :attr:`Universe.atoms` by *index* from the output one would use
   ``u.atoms[index-1]``.


Using the :meth:`HydrogenBondAnalysis.generate_table` method one can reformat
the results as a flat "normalised" table that is easier to import into a
database for further processing. :meth:`HydrogenBondAnalysis.save_table` saves
the table to a pickled file. The table itself is a :class:`numpy.recarray`.


Detection of hydrogen bonds
---------------------------

Hydrogen bonds are recorded based on a geometric criterion:

1. The distance between acceptor and hydrogen is less than or equal to
   *distance* (default is 3 Å).

2. The angle between donor-hydrogen-acceptor is greater than or equal to
   *angle* (default is 120º).

The cut-off values *angle* and *distance* can be set as keywords to
:class:`HydrogenBondAnalysis`.

Donor and acceptor heavy atoms are detected from atom names. The current
defaults are appropriate for the CHARMM27 and GLYCAM06 force fields as defined
in Table `Default atom names for hydrogen bonding analysis`_.

Hydrogen atoms bonded to a donor are searched with one of two algorithms,
selected with the *detect_hydrogens* keyword.

*distance*

   Searches for all hydrogens (name "H*" or name "[123]H" or type "H") in the
   same residue as the donor atom within a cut-off distance of 1.2 Å.

*heuristic*

   Looks at the next three atoms in the list of atoms following the donor and
   selects any atom whose name matches (name "H*" or name "[123]H"). For

The *distance* search is more rigorous but slower and is set as the
default. Until release 0.7.6, only the heuristic search was implemented.

.. versionchanged:: 0.7.6
   Distance search added (see
   :meth:`HydrogenBondAnalysis._get_bonded_hydrogens_dist`) and heuristic
   search improved (:meth:`HydrogenBondAnalysis._get_bonded_hydrogens_list`)

.. _Default atom names for hydrogen bonding analysis:

.. table:: Default heavy atom names for CHARMM27 force field.

   =========== ==============  =========== ====================================
   group       donor           acceptor    comments
   =========== ==============  =========== ====================================
   main chain  N               O
   water       OH2, OW         OH2, OW     SPC, TIP3P, TIP4P (CHARMM27,Gromacs)

   ARG         NE, NH1, NH2
   ASN         ND2             OD1
   ASP                         OD1, OD2
   CYS         SG
   CYH                         SG          possible false positives for CYS
   GLN         NE2             OE1
   GLU                         OE1, OE2
   HIS         ND1, NE2        ND1, NE2    presence of H determines if donor
   HSD         ND1             NE2
   HSE         NE2             ND1
   HSP         ND1, NE2
   LYS         NZ
   MET                         SD          see e.g. [Gregoret1991]_
   SER         OG              OG
   THR         OG1             OG1
   TRP         NE1
   TYR         OH              OH
   =========== ==============  =========== ====================================

.. table:: Heavy atom types for GLYCAM06 force field.

   =========== =========== ==================
   element     donor       acceptor
   =========== =========== ==================
   N           N,NT,N3     N,NT
   O           OH,OW       O,O2,OH,OS,OW,OY
   S                       SM
   =========== =========== ==================

Donor and acceptor names for the CHARMM27 force field will also work for e.g.
OPLS/AA (tested in Gromacs). Residue names in the table are for information
only and are not taken into account when determining acceptors and donors.
This can potentially lead to some ambiguity in the assignment of
donors/acceptors for residues such as histidine or cytosine.

For more information about the naming convention in GLYCAM06 have a look at the
`Carbohydrate Naming Convention in Glycam`.

.. _`Carbohydrate Naming Convention in Glycam`:
   http://glycam.ccrc.uga.edu/documents/FutureNomenclature.htm

The lists of donor and acceptor names can be extended by providing lists of
atom names in the *donors* and *acceptors* keywords to
:class:`HydrogenBondAnalysis`. If the lists are entirely inappropriate
(e.g. when analysing simulations done with a force field that uses very
different atom names) then one should either use the value "other" for *forcefield*
to set no default values, or derive a new class and set the default list oneself::

 class HydrogenBondAnalysis_OtherFF(HydrogenBondAnalysis):
       DEFAULT_DONORS = {"OtherFF": tuple(set([...]))}
       DEFAULT_ACCEPTORS = {"OtherFF": tuple(set([...]))}

Then simply use the new class instead of the parent class and call it with
*forcefield* = "OtherFF". Please also consider to contribute the list of heavy
atom names to MDAnalysis.

.. rubric:: References

.. [Gregoret1991] L.M. Gregoret, S.D. Rader, R.J. Fletterick, and
   F.E. Cohen. Hydrogen bonds involving sulfur atoms in proteins. Proteins,
   9(2):99–107, 1991. `10.1002/prot.340090204`_.

.. _`10.1002/prot.340090204`: http://dx.doi.org/10.1002/prot.340090204


Example
-------

All protein-water hydrogen bonds can be analysed with ::

  import MDAnalysis.analysis.hbonds

  u = MDAnalysis.Universe(PSF, PDB, permissive=True)
  h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'protein', 'resname TIP3', distance=3.0, angle=120.0)
  h.run()

The results are stored as the attribute
:attr:`HydrogenBondAnalysis.timeseries`; see :ref:`Analysis Output`
for the format and further options.

.. Note::

   Due to the way :class:`HydrogenBondAnalysis` is implemented, it is
   more efficient to have the second selection (*selection2*) be the
   *larger* group, e.g. the water when looking at water-protein
   H-bonds or the whole protein when looking at ligand-protein
   interactions.

.. TODO: how to analyse the ouput and notes on selection updating


Classes
-------

.. autoclass:: HydrogenBondAnalysis
   :members:

   .. attribute:: timesteps

      List of the times of each timestep. This can be used together with
      :attr:`~HydrogenBondAnalysis.timeseries` to find the specific time point
      of a hydrogen bond existence, or see :attr:`~HydrogenBondAnalysis.table`.

   .. attribute:: timeseries

      Results of the hydrogen bond analysis, stored for each frame. In
      the following description, # indicates comments that are not
      part of the output::

        results = [
            [ # frame 1
               [ # hbond 1
                  <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
               ],
               [ # hbond 2
                  <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
               ],
               ....
            ],
            [ # frame 2
              [ ... ], [ ... ], ...
            ],
            ...
        ]

      The time of each step is not stored with each hydrogen bond frame but in
      :attr:`~HydrogenBondAnalysis.timesteps`.

      .. Note::

         The *index* is a 1-based index. To get the :attr:`Atom.number` (the
         0-based index typically used in MDAnalysis simply subtract 1. For
         instance, to find an atom in :attr:`Universe.atoms` by *index* one
         would use ``u.atoms[index-1]``.

   .. attribute:: table

      A normalised table of the data in
      :attr:`HydrogenBondAnalysis.timeseries`, generated by
      :meth:`HydrogenBondAnalysis.generate_table`. It is a
      :class:`numpy.recarray` with the following columns:

          0. "time"
          1. "donor_idx"
          2. "acceptor_idx"
          3. "donor_resnm"
          4. "donor_resid"
          5. "donor_atom"
          6. "acceptor_resnm"
          7. "acceptor_resid"
          8. "acceptor_atom"
          9. "distance"
          10. "angle"

      It takes up more space than
      :attr:`~HydrogenBondAnalysis.timeseries` but it is easier to
      analyze and to import into databases (e.g. using recsql_).

      .. Note::

         The *index* is a 1-based index. To get the :attr:`Atom.number` (the
         0-based index typically used in MDAnalysis simply subtract 1. For
         instance, to find an atom in :attr:`Universe.atoms` by *index* one
         would use ``u.atoms[index-1]``.


   .. automethod:: _get_bonded_hydrogens

   .. automethod:: _get_bonded_hydrogens_dist

   .. automethod:: _get_bonded_hydrogens_list

"""

from __future__ import with_statement

from collections import defaultdict
import numpy

from MDAnalysis import MissingDataWarning, NoDataError
from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS
from MDAnalysis.core.util import norm, angle, parse_residue
from MDAnalysis.core.log import ProgressMeter

import warnings
import logging
logger = logging.getLogger('MDAnalysis.analysis.hbonds')

class HydrogenBondAnalysis(object):
    """Perform a hydrogen bond analysis

    The analysis of the trajectory is performed with the
    :meth:`HydrogenBondAnalysis.run` method. The result is stored in
    :attr:`HydrogenBondAnalysis.timeseries`. See
    :meth:`~HydrogenBondAnalysis.run` for the format.

    The default atom names are taken from the CHARMM 27 force field files, which
    will also work for e.g. OPLS/AA in Gromacs, and GLYCAM06.

    *Donors* (associated hydrogens are deduced from topology)
      *CHARMM 27*
        N of the main chain, water OH2/OW, ARG NE/NH1/NH2, ASN ND2, HIS ND1/NE2,
        SER OG, TYR OH, CYS SG, THR OG1, GLN NE2, LYS NZ, TRP NE1
      *GLYCAM06*
        N,NT,N3,OH,OW

    *Acceptors*
      *CHARMM 27*
        O of the main chain, water OH2/OW, ASN OD1, ASP OD1/OD2, CYH SG, GLN OE1,
        GLU OE1/OE2, HIS ND1/NE2, MET SD, SER OG, THR OG1, TYR OH
      *GLYCAM06*
        N,NT,O,O2,OH,OS,OW,OY,P,S,SM

    .. SeeAlso:: Table :ref:`Default atom names for hydrogen bonding analysis`

    .. versionchanged:: 0.7.6
       DEFAULT_DONORS/ACCEPTORS is now embedded in a dict to switch between
       default values for different force fields.
    """

    # use tuple(set()) here so that one can just copy&paste names from the
    # table; set() takes care for removing duplicates. At the end the
    # DEFAULT_DONORS and DEFAULT_ACCEPTORS should simply be tuples.

    #: default heavy atom names whose hydrogens are treated as *donors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`)
    #: Use the keyword *donors* to add a list of additional donor names.
    DEFAULT_DONORS = {'CHARMM27': tuple(set(['N', 'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2',
                                'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH'])),
                      'GLYCAM06': tuple(set(['N','NT','N3','OH','OW'])),
                      'other':  tuple(set([]))}

    #: default atom names that are treated as hydrogen *acceptors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`)
    #: Use the keyword *acceptors* to add a list of additional acceptor names.
    DEFAULT_ACCEPTORS = {'CHARMM27': tuple(set(['O', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1',
                                     'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'])),
                         'GLYCAM06': tuple(set(['N','NT','O','O2','OH','OS','OW','OY','SM'])),
                         'other':  tuple(set([]))}

    #: A :class:`collections.defaultdict` of covalent radii of common donors
    #: (used in :meth`_get_bonded_hydrogens_list` to check if a hydrogen is
    #: sufficiently close to its donor heavy atom). Values are stored for
    #: N, O, P, and S. Any other heavy atoms are assumed to have hydrogens
    #: covalently bound at a maximum distance of 1.5 Å.
    r_cov = defaultdict(lambda : 1.5,   # default value
                        N=1.31, O=1.31, P=1.58, S=1.55)

    def __init__(self, universe, selection1='protein', selection2='all', selection1_type='both',
                 update_selection1=False, update_selection2=False, filter_first=True,
                 distance=3.0, angle=120.0,
                 forcefield='CHARMM27', donors=None, acceptors=None,
                 start=None, stop=None, step=None,
                 verbose=False, detect_hydrogens='distance'):
        """Set up calculation of hydrogen bonds between two selections in a universe.

        :Arguments:
          *universe*
            Universe object
          *selection1*
            Selection string for first selection ['protein']
          *selection2*
            Selection string for second selection ['all']
          *selection1_type*
            Selection 1 can be 'donor', 'acceptor' or 'both' ['both']
          *update_selection1*
            Update selection 1 at each frame? [``False``]
          *update_selection2*
            Update selection 2 at each frame? [``False``]
          *filter_first*
            Filter selection 2 first to only atoms 3*distance away [``True``]
          *distance*
            Distance cutoff for hydrogen bonds; only interactions with a H-A distance
            <= *distance* (and the appropriate D-H-A angle, see *angle*) are
            recorded. [3.0
          *angle*
            Angle cutoff for hydrogen bonds; an ideal H-bond has an angle of
            180º.  A hydrogen bond is only recorded if the D-H-A angle is
            >=  *angle*. The default of 120º also finds fairly non-specific
            hydrogen interactions and a possibly better value is 150º. [120.0]
          *forcefield*
            Name of the forcefield used. Switches between different
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS` and
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS` values.
            Available values: "CHARMM27", "GLYCAM06", "other" ["CHARMM27"]
          *donors*
            Extra H donor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS`), must be a sequence.
          *acceptors*
            Extra H acceptor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS`), must be a sequence.
          *start*
            starting frame-index for analysis, ``None`` is the first one, 0.
            *start* and *stop* are 0-based frame indices and are used to slice
            the trajectory (if supported) [``None``]
          *stop*
            last trajectory frame for analysis, ``None`` is the last one [``None``]
          *step*
            read every *step* between *start* and *stop*, ``None`` selects 1.
            Note that not all trajectory readers perform well with a step different
            from 1 [``None``]
          *verbose*
            If set to ``True`` enables per-frame debug logging. This is disabled
            by default because it generates a very large amount of output in
            the log file. (Note that a logger must have been started to see
            the output, e.g. using :func:`MDAnalysis.start_logging`.)
          *detect_hydrogens*
            Determine the algorithm to find hydrogens connected to donor
            atoms. Can be "distance" (default; finds all hydrogens in the
            donor's residue within a cutoff of the donor) or "heuristic"
            (looks for the next few atoms in the atom list). "distance" should
            always give the correct answer but "heuristic" is faster,
            especially when the donor list is updated each
            for each frame. ["distance"]

        The timeseries is accessible as the attribute :attr:`HydrogenBondAnalysis.timeseries`.

        .. versionchanged:: 0.7.6
           New *verbose* keyword (and per-frame debug logging disabled by
           default).

           New *detect_hydrogens* keyword to switch between two different
           algorithms to detect hydrogens bonded to donor. "distance" is a new,
           rigorous distance search within the residue of the donor atom,
           "heuristic" is the previous list scan (improved with an additional
           distance check).

           New *forcefield* keyword to switch between different values of
           DEFAULT_DONORS/ACCEPTORS to accomodate different force fields.
           Also has an option "other" for no default values.
        """
        self._get_bonded_hydrogens_algorithms = {
            "distance": self._get_bonded_hydrogens_dist,      # 0.7.6 default
            "heuristic": self._get_bonded_hydrogens_list,     # pre 0.7.6
            }
        if not detect_hydrogens in self._get_bonded_hydrogens_algorithms:
            raise ValueError("detect_hydrogens must be one of %r" %
                             self._get_bonded_hydrogens_algorithms.keys())
        self.detect_hydrogens = detect_hydrogens

        self.u = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection1_type = selection1_type
        self.update_selection1 = update_selection1
        self.update_selection2 = update_selection2
        self.filter_first = filter_first
        self.distance = distance
        self.angle = angle
        self.traj_slice = slice(start if isinstance(start, int) else None, # internal frames are 0 based
                                stop if isinstance(stop, int) else None,
                                step)

        # set up the donors/acceptors lists
        if donors is None:
            donors = []
        if acceptors is None:
            acceptors = []
        self.donors = tuple(set(self.DEFAULT_DONORS[forcefield]).union(donors))
        self.acceptors = tuple(set(self.DEFAULT_ACCEPTORS[forcefield]).union(acceptors))

        if not (self.selection1 and self.selection2):
            raise ValueError('HydrogenBondAnalysis: invalid selections')
        elif self.selection1_type not in ('both', 'donor', 'acceptor'):
            raise ValueError('HydrogenBondAnalysis: Invalid selection type %s' % self.selection1_type)

        self.timeseries = None  # final result
        self.timesteps = None  # time for each frame

        self.table = None # placeholder for output table

        self.verbose = True # always enable debug output for initial selection update
        self._update_selection_1()
        self._update_selection_2()
        self.verbose = verbose  # per-frame debugging output?

    def _get_bonded_hydrogens(self, atom, **kwargs):
        """Find hydrogens bonded to *atom*.

        This method is typically not called by a user but it is documented to
        facilitate understanding of the internals of
        :class:`HydrogenBondAnalysis`.

        :Returns: list of hydrogens (can be a
                  :class:`~MDAnalysis.core.AtomGroup.AtomGroup`) or empty list
                  ``[]`` if none were found.

        .. SeeAlso::

           :meth:`_get_bonded_hydrogens_dist` and :meth:`_get_bonded_hydrogens_list`

        .. versionchanged:: 0.7.6
           Can switch algorithm by using the *detect_hydrogens* keyword to the
           constructor. *kwargs* can be used to supply arguments for algorithm.
        """
        return self._get_bonded_hydrogens_algorithms[self.detect_hydrogens](atom, **kwargs)

    def _get_bonded_hydrogens_dist(self, atom):
        """Find hydrogens bonded within *cutoff* to *atom*.

        * hydrogens are detected by either name ("H*", "[123]H*") or type
          ("H"); this is not fool-proof as the atom type is not always a
          character but the name pattern should catch most typical occurrences.

        * The distance from *atom* is calculated for all hydrogens in the
          residue and only those within a cutoff are kept. The cutoff depends
          on the heavy atom (more precisely, on its element, which is taken as
          the first letter of its name ``atom.name[0]``) and is parameterized
          in :attr:`HydrogenBondAnalysis.r_cov`. If no match is found then the
          default of 1.5 Å is used.

        The performance of this implementation could be improved once the
        topology always contains bonded information; it currently uses the
        selection parser with an "around" selection.

        .. versionadded:: 0.7.6
        """
        try:
            return atom.residue.selectAtoms("(name H* or name 1H* or name 2H* or name 3H* or type H) and around %f name %s" %
                                            (self.r_cov[atom.name[0]], atom.name))
        except NoDataError:
            return []

    def _get_bonded_hydrogens_list(self, atom, **kwargs):
        """Find "bonded" hydrogens to the donor *atom*.

        At the moment this relies on the **assumption** that the
        hydrogens are listed directly after the heavy atom in the
        topology. If this is not the case then this function will
        fail.

        Hydrogens are detected by name ``H*``, ``[123]H*`` and they have to be
        within a maximum distance from the heavy atom. The cutoff distance
        depends on the heavy atom and is parameterized in
        :attr:`HydrogenBondAnalysis.r_cov`.

        .. versionchanged:: 0.7.6

           Added detection of ``[123]H`` and additional check that a
           selected hydrogen is bonded to the donor atom (i.e. its
           distance to the donor is less than the covalent radius
           stored in :attr:`HydrogenBondAnalysis.r_cov` or the default
           1.5 Å).

           Changed name to
           :meth:`~HydrogenBondAnalysis._get_bonded_hydrogens_list`
           and added *kwargs* so that it can be used instead of
           :meth:`~HydrogenBondAnalysis._get_bonded_hydrogens_dist`.

        """
        warnings.warn("_get_bonded_hydrogens_list() does not always find "
                      "all hydrogens; detect_hydrogens='distance' is safer.",
                      category=DeprecationWarning)
        try:
            hydrogens = [a for a in self.u.atoms[atom.number+1:atom.number+4]
                         if a.name.startswith(('H','1H','2H','3H')) \
                             and self.calc_eucl_distance(atom,a) < self.r_cov[atom.name[0]]]
        except IndexError:
            hydrogens = []  # weird corner case that atom is the last one in universe
        return hydrogens

    def _update_selection_1(self):
        self._s1 = self.u.selectAtoms(self.selection1)
        self.logger_debug("Size of selection 1: %d atoms" % len(self._s1))
        self._s1_donors = {}
        self._s1_donors_h = {}
        self._s1_acceptors = {}
        if self.selection1_type in ('donor', 'both'):
            self._s1_donors = self._s1.selectAtoms(' or '.join([ 'name %s' % i for i in self.donors ]))
            self._s1_donors_h = {}
            for i,d in enumerate(self._s1_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._s1_donors_h[i] = tmp
            self.logger_debug("Selection 1 donors: %d" % len(self._s1_donors))
            self.logger_debug("Selection 1 donor hydrogens: %d" % len(self._s1_donors_h))
        if self.selection1_type in ('acceptor', 'both'):
            self._s1_acceptors = self._s1.selectAtoms(' or '.join([ 'name %s' % i for i in self.acceptors ]))
            self.logger_debug("Selection 1 acceptors: %d" % len(self._s1_acceptors))

    def _update_selection_2(self):
        self._s2 = self.u.selectAtoms(self.selection2)
        if self.filter_first:
            self.logger_debug("Size of selection 2 before filtering: %d atoms" % len(self._s2))
            ns_selection_2 = NS.AtomNeighborSearch(self._s2)
            self._s2 = ns_selection_2.search_list(self._s1, 3.*self.distance)
        self.logger_debug("Size of selection 2: %d atoms" % len(self._s2))

        self._s2_donors = {}
        self._s2_donors_h = {}
        self._s2_acceptors = {}
        if self.selection1_type in ('donor', 'both'):
            self._s2_acceptors = self._s2.selectAtoms(' or '.join([ 'name %s' % i for i in self.acceptors ]))
            self.logger_debug("Selection 2 acceptors: %d" % len(self._s2_acceptors))
        if self.selection1_type in ('acceptor', 'both'):
            self._s2_donors = self._s2.selectAtoms(' or '.join([ 'name %s' % i for i in self.donors ]))
            self._s2_donors_h = {}
            for i,d in enumerate(self._s2_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._s2_donors_h[i] = tmp
            self.logger_debug("Selection 2 donors: %d" % len(self._s2_donors))
            self.logger_debug("Selection 2 donor hydrogens: %d" % len(self._s2_donors_h))

    def logger_debug(self, *args):
        if self.verbose:
            logger.debug(*args)

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries.

        Stores the hydrogen bond data per frame
        as :attr:`HydrogenBondAnalysis.timeseries` (see there for
        output format).

        .. SeeAlso:: :meth:`HydrogenBondAnalysis.generate_table` for processing
                     the data into a different format.

        .. versionchanged:: 0.7.6
           Results are not returned, only stored in
           :attr:`~HydrogenBondAnalysis.timeseries` and duplicate hydrogen bonds
           are removed from output (can be suppressed with *remove_duplicates* =
           ``False``)

        """
        logger.info("HBond analysis: starting")
        logger.debug("HBond analysis: donors    %r", self.donors)
        logger.debug("HBond analysis: acceptors %r", self.acceptors)

        remove_duplicates = kwargs.pop('remove_duplicates', True) # False: old behaviour
        if not remove_duplicates:
            logger.warn("Hidden feature remove_duplicates = True activated: you will probably get duplicate H-bonds.")

        verbose = kwargs.pop('verbose', False)
        if verbose != self.verbose:
            self.verbose = verbose
            logger.debug("Toggling verbose to %r", self.verbose)
        if not self.verbose:
            logger.debug("HBond analysis: For full step-by-step debugging output use verbose=True")

        self.timeseries = []
        self.timesteps = []

        logger.info("checking trajectory...")  # numframes can take a while!
        try:
            frames = numpy.arange(self.u.trajectory.numframes)[self.traj_slice]
        except:
            logger.error("Problem reading trajectory or trajectory slice incompatible.")
            logger.exception()
            raise
        pm = ProgressMeter(len(frames),
                           format="HBonds frame %(step)5d/%(numsteps)d [%(percentage)5.1f%%]\r")

        try:
            self.u.trajectory.time
            def _get_timestep():
                return self.u.trajectory.time
            logger.debug("HBond analysis is recording time step")
        except NotImplementedError:
            # chained reader or xyz(?) cannot do time yet
            def _get_timestep():
                return self.u.trajectory.frame
            logger.warn("HBond analysis is recording frame number instead of time step")


        logger.info("Starting analysis (frame index start=%d stop=%d, step=%d)",
                    (self.traj_slice.start or 0),
                    (self.traj_slice.stop or self.u.trajectory.numframes), self.traj_slice.step or 1)

        for ts in self.u.trajectory[self.traj_slice]:
            # all bonds for this timestep
            frame_results = []
            # dict of tuples (atomid, atomid) for quick check if
            # we already have the bond (to avoid duplicates)
            already_found = {}

            frame = ts.frame
            timestep = _get_timestep()
            self.timesteps.append(timestep)

            pm.echo(ts.frame)
            self.logger_debug("Analyzing frame %(frame)d, timestep %(timestep)f ps", vars())
            if self.update_selection1:
                self._update_selection_1()
            if self.update_selection2:
                self._update_selection_2()

            if self.selection1_type in ('donor', 'both'):
                self.logger_debug("Selection 1 Donors <-> Acceptors")
                ns_acceptors = NS.AtomNeighborSearch(self._s2_acceptors)
                for i,donor_h_set in self._s1_donors_h.items():
                    d = self._s1_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search_list(AtomGroup([h]), self.distance)
                        for a in res:
                            angle = self.calc_angle(d,h,a)
                            dist = self.calc_eucl_distance(h,a)
                            if angle >= self.angle and dist <= self.distance:
                                self.logger_debug("S1-D: %s <-> S2-A: %s %f A, %f DEG" % (h.number+1, a.number+1, dist, angle))
                                #self.logger_debug("S1-D: %r <-> S2-A: %r %f A, %f DEG" % (h, a, dist, angle))
                                frame_results.append([h.number+1, a.number+1, '%s%s:%s' % (h.resname, repr(h.resid), h.name), '%s%s:%s' % (a.resname, repr(a.resid), a.name), dist, angle])
                                already_found[(h.number+1, a.number+1)] = True
            if self.selection1_type in ('acceptor', 'both'):
                self.logger_debug("Selection 1 Acceptors <-> Donors")
                ns_acceptors = NS.AtomNeighborSearch(self._s1_acceptors)
                for i,donor_h_set in self._s2_donors_h.items():
                    d = self._s2_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search_list(AtomGroup([h]), self.distance)
                        for a in res:
                            if remove_duplicates and \
                                    ((h.number+1, a.number+1) in already_found or \
                                         (a.number+1, h.number+1) in already_found):
                                continue
                            angle = self.calc_angle(d,h,a)
                            dist = self.calc_eucl_distance(h,a)
                            if angle >= self.angle and dist <= self.distance:
                                self.logger_debug("S1-A: %s <-> S2-D: %s %f A, %f DEG" % (a.number+1, h.number+1, dist, angle))
                                #self.logger_debug("S1-A: %r <-> S2-D: %r %f A, %f DEG" % (a, h, dist, angle))
                                frame_results.append([h.number+1, a.number+1, '%s%s:%s' % (h.resname, repr(h.resid), h.name), '%s%s:%s' % (a.resname, repr(a.resid), a.name), dist, angle])
            self.timeseries.append(frame_results)

        logger.info("HBond analysis: complete; timeseries with %d hbonds in %s.timeseries",
                    self.count_by_time().count.sum(), self.__class__.__name__)

    def calc_angle(self, d, h, a):
        """Calculate the angle (in degrees) between two atoms with H at apex."""
        v1 = h.pos-d.pos
        v2 = h.pos-a.pos
        if numpy.all(v1 == v2):
            return 0.0
        return numpy.rad2deg(angle(v1, v2))

    def calc_eucl_distance(self, a1, a2):
        """Calculate the Euclidean distance between two atoms. """
        return norm(a2.pos - a1.pos)

    def generate_table(self):
        """Generate a normalised table of the results.

        The table is stored as a :class:`numpy.recarray` in the
        attribute :attr:`~HydrogenBondAnalysis.table` and can be used
        with e.g. `recsql`_.

        Columns:
          0. "time"
          1. "donor_idx"
          2. "acceptor_idx"
          3. "donor_resnm"
          4. "donor_resid"
          5. "donor_atom"
          6. "acceptor_resnm"
          7. "acceptor_resid"
          8. "acceptor_atom"
          9. "distance"
          10. "angle"

        .. _recsql: http://pypi.python.org/pypi/RecSQL
        """
        from itertools import izip
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        num_records = numpy.sum([len(hframe) for hframe in self.timeseries])
        dtype = [("time",float), ("donor_idx",int), ("acceptor_idx",int),
                 ("donor_resnm","|S4"), ("donor_resid",int), ("donor_atom","|S4"),
                 ("acceptor_resnm","|S4"), ("acceptor_resid",int), ("acceptor_atom","|S4"),
                 ("distance",float), ("angle",float)]
        self.table = numpy.recarray((num_records,), dtype=dtype)

        # according to Lukas' notes below, using a recarray at this stage is ineffective
        # and speedups of ~x10 could be achieved by filling a standard array
        # (perhaps at the cost of less clarity... but that might just be my code ;-) -- orbeckst)
        cursor = 0  # current row
        for t,hframe in izip(self.timesteps, self.timeseries):
            if len(hframe) == 0:
                continue    # not really necessary, should also work without
            self.table[cursor:cursor+len(hframe)].time = t
            for donor_idx, acceptor_idx, donor, acceptor, distance, angle in hframe:
                r = self.table[cursor]
                r.donor_idx = donor_idx
                r.donor_resnm, r.donor_resid, r.donor_atom = parse_residue(donor)
                r.acceptor_idx = acceptor_idx
                r.acceptor_resnm, r.acceptor_resid, r.acceptor_atom = parse_residue(acceptor)
                r.distance = distance
                r.angle = angle
                cursor += 1
        assert cursor == num_records, "Internal Error: Not all HB records stored"
        logger.debug("HBond: Stored results as table with %(num_records)d entries.", vars())

    def save_table(self, filename="hbond_table.pickle"):
        """Saves :attr:`~HydrogenBondAnalysis.table` to a pickled file.

        Load with ::

           import cPickle
           table = cPickle.load(open(filename))

        .. SeeAlso:: :mod:`cPickle` module and :class:`numpy.recarray`
        """
        import cPickle
        if self.table is None:
            self.generate_table()
        cPickle.dump(self.table, open(filename, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL)

    def count_by_time(self):
        """Counts the number of hydrogen bonds per timestep.

        :Returns: a class:`numpy.recarray`
        """
        from itertools import izip, imap
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        out = numpy.empty((len(self.timesteps),), dtype=[('time', float), ('count', int)])
        for cursor,time_count in enumerate(izip(self.timesteps, imap(len, self.timeseries))):
            out[cursor] = time_count
        return out.view(numpy.recarray)

    def count_by_type(self):
        """Counts the frequency of hydrogen bonds of a specific type.

        Processes :attr:`HydrogenBondAnalysis.timeseries` and returns
        a :class:`numpy.recarray` containing atom indices, residue
        names, residue numbers (for donors and acceptors) and the
        fraction of the total time during which the hydrogen bond was
        detected.

        :Returns: a class:`numpy.recarray`
        """
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        hbonds = defaultdict(int)
        for hframe in self.timeseries:
            for donor_idx, acceptor_idx, donor, acceptor, distance, angle in hframe:
                donor_resnm, donor_resid, donor_atom = parse_residue(donor)
                acceptor_resnm, acceptor_resid, acceptor_atom = parse_residue(acceptor)
                # generate unambigous key for current hbond
                # (the donor_heavy_atom placeholder '?' is added later)
                hb_key = (donor_idx, acceptor_idx,
                          donor_resnm,    donor_resid,  "?", donor_atom,
                          acceptor_resnm, acceptor_resid, acceptor_atom)

                hbonds[hb_key] += 1

        # build empty output table
        dtype = [('donor_idx', int), ('acceptor_idx', int),
                ('donor_resnm', 'S4'), ('donor_resid', int), ('donor_heavy_atom', 'S4'), ('donor_atom', 'S4'),
                ('acceptor_resnm', 'S4'), ('acceptor_resid', int), ('acceptor_atom', 'S4'),
                ('frequency', float)]
        out = numpy.empty((len(hbonds),), dtype=dtype)

        # float because of division later
        tsteps = float(len(self.timesteps))
        for cursor, (key, count) in enumerate(hbonds.iteritems()):
            out[cursor] = key + (count/tsteps,)

        # return array as recarray
        # The recarray has not been used within the function, because accessing the
        # the elements of a recarray (3.65 us) is much slower then accessing those
        # of a ndarray (287 ns).
        r = out.view(numpy.recarray)

        # patch in donor heavy atom names (replaces '?' in the key)
        h2donor = self._donor_lookup_table_byindex()
        r.donor_heavy_atom[:] = [h2donor[idx-1] for idx in r.donor_idx]

        return r

    def timesteps_by_type(self):
        """Frames during which each hydrogen bond existed, sorted by hydrogen bond.

        Processes :attr:`HydrogenBondAnalysis.timeseries` and returns
        a :class:`numpy.recarray` containing atom indices, residue
        names, residue numbers (for donors and acceptors) and a list
        of timesteps at which the hydrogen bond was detected.

        :Returns: a class:`numpy.recarray`
        """
        from itertools import izip
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        hbonds = defaultdict(list)
        for (t,hframe) in izip(self.timesteps, self.timeseries):
            for donor_idx, acceptor_idx, donor, acceptor, distance, angle in hframe:
                donor_resnm, donor_resid, donor_atom = parse_residue(donor)
                acceptor_resnm, acceptor_resid, acceptor_atom = parse_residue(acceptor)
                # generate unambigous key for current hbond
                # (the donor_heavy_atom placeholder '?' is added later)
                hb_key = (donor_idx, acceptor_idx,
                          donor_resnm,    donor_resid,  "?", donor_atom,
                          acceptor_resnm, acceptor_resid, acceptor_atom)
                hbonds[hb_key].append(t)

        out_nrows = 0
        # count number of timesteps per key to get length of output table
        for ts_list in hbonds.itervalues():
            out_nrows += len(ts_list)

        # build empty output table
        dtype = [('donor_idx', int), ('acceptor_idx', int),
                ('donor_resnm', 'S4'),    ('donor_resid', int), ('donor_heavy_atom', 'S4'),  ('donor_atom', 'S4'),
                ('acceptor_resnm', 'S4'), ('acceptor_resid', int), ('acceptor_atom', 'S4'),
                ('time', float)]
        out = numpy.empty((out_nrows,), dtype=dtype)

        out_row = 0
        for (key, times) in hbonds.iteritems():
            for tstep in times:
                out[out_row] = key + (tstep,)
                out_row += 1

        # return array as recarray
        # The recarray has not been used within the function, because accessing the
        # the elements of a recarray (3.65 us) is much slower then accessing those
        # of a ndarray (287 ns).
        r = out.view(numpy.recarray)

        # patch in donor heavy atom names (replaces '?' in the key)
        h2donor = self._donor_lookup_table_byindex()
        r.donor_heavy_atom[:] = [h2donor[idx-1] for idx in r.donor_idx]

        return r

    def _donor_lookup_table_byres(self):
        """Look-up table to identify the donor heavy atom from resid and hydrogen name.

        Assumptions:
        * resids are unique
        * hydrogen atom names are unique within a residue
        * selections have not changed (because we are simply looking at the last content
          of the donors and donor hydrogen lists)

        Donors from *selection1* and *selection2* are merged.

        Output dictionary ``h2donor`` can be used as::

           heavy_atom_name = h2donor[resid][hydrogen_name]

        """
        s1d = self._s1_donors    # list of donor Atom instances
        s1h = self._s1_donors_h  # dict indexed by donor position in donor list, containg AtomGroups of H
        s2d = self._s2_donors
        s2h = self._s2_donors_h

        def _make_dict(donors, hydrogens):
            # two steps so that entry for one residue can be UPDATED for multiple donors
            d = dict((donors[k].resid, {})  for k in xrange(len(donors)) if k in hydrogens)
            for k in xrange(len(donors)):
                if k in hydrogens:
                    d[donors[k].resid].update(dict((atom.name, donors[k].name) for atom in hydrogens[k]))
            return d

        h2donor = _make_dict(s2d, s2h)    # 2 is typically the larger group
        # merge (in principle h2donor.update(_make_dict(s1d, s1h) should be sufficient
        # with our assumptions but the following should be really safe)
        for resid,names in _make_dict(s1d, s1h).items():
            if resid in h2donor:
                h2donor[resid].update(names)
            else:
                h2donor[resid] = names

        return h2donor

    def _donor_lookup_table_byindex(self):
        """Look-up table to identify the donor heavy atom from hydrogen atom index.

        Assumptions:
        * selections have not changed (because we are simply looking at the last content
          of the donors and donor hydrogen lists)

        Donors from *selection1* and *selection2* are merged.

        Output dictionary ``h2donor`` can be used as::

           heavy_atom_name = h2donor[index]

        .. Note::

           *index* is the 0-based MDAnalysis index
           (:attr:`MDAnalysis.core.AtomGroup.Atom.number`).  The
           tables generated by :class:`HydrogenBondAnalysis` contain
           1-based indices.

        """
        s1d = self._s1_donors    # list of donor Atom instances
        s1h = self._s1_donors_h  # dict indexed by donor position in donor list, containg AtomGroups of H
        s2d = self._s2_donors
        s2h = self._s2_donors_h

        def _make_dict(donors, hydrogens):
            #return dict(flatten_1([(atom.id, donors[k].name) for atom in hydrogens[k]] for k in xrange(len(donors)) if k in hydrogens))
            x = []
            for k in xrange(len(donors)):
                if k in hydrogens:
                    x.extend([(atom.number, donors[k].name) for atom in hydrogens[k]])
            return dict(x)

        h2donor = _make_dict(s2d, s2h)    # 2 is typically the larger group
        h2donor.update(_make_dict(s1d, s1h))

        return h2donor


