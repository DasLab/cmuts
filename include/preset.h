/* preset.h -- the minimap2 presets that reads are aligned under.
 *
 * Two places need the presets: the command line offers them, and the alignment run passes
 * one to minimap2. Both share the cli_choice table below, which is why this header
 * includes cli.h where no other domain header does.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>

#include "cli.h"

/* Every preset that aligns sequencing reads against a reference, which is the alignment
 * cmuts hmm counts. */
typedef enum {
    PRESET_SR,
    PRESET_MAP_ONT,
    PRESET_MAP_HIFI,
    PRESET_MAP_PB,
    PRESET_MAP_ICLR,
    PRESET_LR_HQ,
    /* No preset was given. The option is required, so a run never sees this value. */
    PRESET_UNSET = -1,
} preset_id;

/* The presets that the command line accepts. The list ends with a NULL name. */
extern const cli_choice PRESET_CHOICES[];

/* Returns the name that minimap2 accepts for -x, or NULL for PRESET_UNSET. */
const char *preset_name(preset_id id);

/* Returns true if the sequencing technology produces paired reads. */
bool preset_accepts_pairs(preset_id id);
