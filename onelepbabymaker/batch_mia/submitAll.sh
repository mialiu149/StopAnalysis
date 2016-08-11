#!/bin/bash
voms-proxy-init -voms cms -valid 240:00
condor_submit configs_V80_signalscanv4/condor_V80_signalscanv4_SMS_tchwh_lnbb.cmd
