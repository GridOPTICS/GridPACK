Python API Reference
====================

This chapter documents the module behind the ``gridpack`` command and the
pybind11 bindings added for contingency analysis.

.. module:: gridpack.cli.main

CLI Module
----------

``gridpack.cli.main`` is a thin argparse shell over the high-level
wrappers (``gridpack.PowerFlow``, ``gridpack.DynamicSim``,
``gridpack.ContingencyAnalysis`` and so on); it holds no solver logic of
its own.  The ``gridpack`` console script points at :func:`main`.

Entry Point
~~~~~~~~~~~

.. function:: main(argv=None)

   Parse the command line, open one ``gridpack.Session``, dispatch to the
   subcommand handler and close the session.

   :param argv: Command-line arguments. Defaults to ``sys.argv[1:]``.
   :type argv: list[str] or None
   :returns: Exit code.
   :rtype: int

   Exit codes:

   - ``0`` success
   - ``1`` unexpected error (set ``GRIDPACK_TRACEBACK=1`` to keep the trace)
   - ``2`` bad arguments or unreadable config
   - ``3`` the analysis did not converge

Subcommand Handlers
~~~~~~~~~~~~~~~~~~~

Each subcommand is a function taking the parsed ``argparse.Namespace`` and
the open session.

.. function:: cmd_powerflow(args, session)

   AC power flow.  ``--solver`` overrides the XML ``UseNonLinear``;
   ``--export-psse FORMAT FILE`` writes PSS/E v23, v33 or v34.

.. function:: cmd_dsf(args, session)

   Dynamic simulation, the ``dsf2.x`` path.  ``--output`` is refused: the
   run writes through ``printf`` that no redirect captures.

.. function:: cmd_se(args, session)

   State estimation.  ``--output`` is refused for the same reason.

.. function:: cmd_hadrec(args, session)

   HADREC remedial action control.

.. function:: cmd_emt(args, session)

   Electromagnetic transient simulation.

.. function:: cmd_ca(args, session)

   N-1 contingency analysis.  ``--vlimits VMIN VMAX`` overrides the XML
   voltage band; ``--no-print-calcs`` suppresses per-contingency files.


Contingency Analysis Bindings
-----------------------------

The following classes, methods, and constants were added to the
``gridpack.powerflow`` module to support contingency analysis
through the Python interface.

Contingency Class
~~~~~~~~~~~~~~~~~

.. class:: gridpack.powerflow.Contingency

   Data container describing a single contingency event. Instances
   of this class are passed to
   :meth:`~gridpack.powerflow.Powerflow.setContingency` and
   :meth:`~gridpack.powerflow.Powerflow.unSetContingency`.

   .. attribute:: p_type

      Contingency type. Use the module-level constants
      ``gridpack.powerflow.Generator`` (0) or
      ``gridpack.powerflow.Branch`` (1).

      :type: int

   .. attribute:: p_name

      Human-readable contingency identifier.

      :type: str

   **Line contingency attributes:**

   .. attribute:: p_from

      List of from-bus numbers for lines to be tripped.

      :type: list[int]

   .. attribute:: p_to

      List of to-bus numbers for lines to be tripped.

      :type: list[int]

   .. attribute:: p_ckt

      List of circuit identifiers for lines to be tripped
      (two-character strings, e.g., ``"BL"``).

      :type: list[str]

   .. attribute:: p_saveLineStatus

      Line status flags saved prior to contingency application.
      Populated automatically by ``setContingency()``.

      :type: list[bool]

   **Generator contingency attributes:**

   .. attribute:: p_busid

      List of bus numbers for generators to be tripped.

      :type: list[int]

   .. attribute:: p_genid

      List of generator identifiers to be tripped
      (two-character strings, e.g., ``"1 "``).

      :type: list[str]

   .. attribute:: p_saveGenStatus

      Generator status flags saved prior to contingency application.
      Populated automatically by ``setContingency()``.

      :type: list[bool]

Powerflow Methods
~~~~~~~~~~~~~~~~~

The following methods were added to
:class:`gridpack.powerflow.Powerflow` to support the contingency
analysis workflow:

.. method:: gridpack.powerflow.Powerflow.setContingency(event)

   Apply a contingency to the network by modifying the online status
   of the specified lines or generators. The pre-contingency status
   is saved in the ``event`` object for later restoration.

   :param event: Contingency definition.
   :type event: :class:`gridpack.powerflow.Contingency`
   :returns: ``True`` if the contingency elements were found in the network.
   :rtype: bool

.. method:: gridpack.powerflow.Powerflow.unSetContingency(event)

   Restore the network to its pre-contingency state using the saved
   status flags in the ``event`` object.

   :param event: Contingency to reverse.
   :type event: :class:`gridpack.powerflow.Contingency`
   :returns: ``True`` if the restoration was successful.
   :rtype: bool

.. method:: gridpack.powerflow.Powerflow.writeCABus()

   Write bus-level results formatted for contingency analysis reports.

.. method:: gridpack.powerflow.Powerflow.writeCABranch()

   Write branch-level results formatted for contingency analysis
   reports.

.. method:: gridpack.powerflow.Powerflow.checkSlackCapacity()

   Verify that the slack bus generator's real power output does not
   exceed its rated capacity.

   :returns: ``True`` if the slack generator is within its capacity limits.
   :rtype: bool

.. method:: gridpack.powerflow.Powerflow.getIslandCount()

   Determine the number of electrically connected islands in the
   network. A value greater than 1 indicates that a contingency
   has caused islanding.

   :returns: Number of islands (1 for a fully connected network).
   :rtype: int

.. method:: gridpack.powerflow.Powerflow.hasLoneBus()

   Check whether any buses in the network have become completely
   isolated (no active branch connections).

   :returns: ``True`` if one or more isolated buses exist.
   :rtype: bool

Constants
~~~~~~~~~

.. data:: gridpack.powerflow.Branch

   Integer constant identifying a branch (line) contingency type.
   Value: ``1``.

.. data:: gridpack.powerflow.Generator

   Integer constant identifying a generator contingency type.
   Value: ``0``.
