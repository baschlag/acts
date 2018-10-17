// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

/// @brief Actor as configurator of the Stepper for working with material interaction. It sets up a tracking of the material and the particles mass and adds further configuration properties. 
/// The additional properties will only be updated (fully) if the corresponding flag is set by the user.
struct StepActor
{ 
	// Configurations for Stepper
    /// Boolean flag for energy loss while stepping
    bool m_energyLossFlag = true;
    /// Tolerance for the error of the integration
	double m_tolerance = 5e-5;
    /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
	bool m_includeGgradient = false;
	/// Cut-off value for the momentum
    double m_momentumCutOff = 0.;
    /// User defined flag for updating parameters in the stepper
    bool m_update = false;
    
    /// @brief Main call operator for setting up stepper properties
    ///
    /// @tparam propagator_state_t Type of the propagator state
    ///
    /// @param [in, out] state State of the propagator
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state) const
  {
    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

	// Initialize all parameters
	if(state.stepping.pathAccumulated == 0.)
	{
		// Let the stepper track the volume and particles mass
		state.stepping.volume = &state.navigation.currentVolume;
		state.stepping.mass = &(state.options.mass);
		
		// Initialize user defined parameters
		state.stepping.energyLossFlag = m_energyLossFlag;
		state.stepping.tolerance = m_tolerance;
		state.stepping.includeGgradient = m_includeGgradient;
		state.stepping.momentumCutOff = m_momentumCutOff;
	}
	else
	{
		// Walk over every updateable parameter and perform the update if wished 
		// (but also as often as wished)
		if(m_update)
		{
			state.stepping.energyLossFlag = m_energyLossFlag;
			state.stepping.tolerance = m_tolerance;
			state.stepping.includeGgradient = m_includeGgradient;
			state.stepping.momentumCutOff = m_momentumCutOff;
		}
	}
  }  
};
}
