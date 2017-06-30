Developer Notes / Tools
=======================

Assorted "How-To's" for developers

How to contribute changes
-------------------------
- Clone the repository if you have write access to the main repo, fork the repository if you are a collaborator.
- Make a new branch with `git checkout -b {your branch name}`
- Make changes and test your code
- Push the branch to the repo (either the main or your fork) with `git push -u origin {your branch name}`
  * Note that `origin` is the default name assigned to the remote, yours may be different
- Make a PR on GitHub with your changes
- We'll review the changes and get your code into the repo after lively discussion!


Checklist for all Updates
-------------------------
- [ ] Update `setup.py` version number (see specific update type for details)
- [ ] Make sure there is an/are issue(s) opened for your specific update
- [ ] Create the PR, referencing the issue
- [ ] Debug the PR as needed until tests pass
- [ ] Tag the final, debugged version as the one in `setup.py`
   *  `git tag -a X.Y.Z [latest pushed commit] && git push --follow-tags`
- [ ] Get the PR merged in


Checklist for Major Revisions (`openmm-forcefields` X.Y+1.0)
--------------------------------------------------------------
- [ ] Make sure all issues related to the milestone will be closed by this commit or moved to future releases
- [ ] Update `docs/whatsnew.rst`
- [ ] Update `setup.py` with version number and `ISRELEASED` to `True`
- [ ] Do the steps for All Upates
- [ ] Create a new release on GitHub, reference the tag and copy the changes in `docs/whatsnew.rst`
- [ ] Update the `omnia-md/conda-recipies` repo by creating a new PR with updated versions
  * Be sure to pin dependencies to fixed version

Checklist for Minor Revisions (`openmm-forcefields` X.Y.Z+1)
--------------------------------------------
- [ ] Update `setup.py` with the correct Z version number in X.Y.Z
- [ ] In `setup.py`, set `ISRELEASED` to `False`
- [ ] Do all the steps for All Updates
  * If this is a critical bugfix (i.e. `openmm-forcefields` X.Y.0 is broken and/or flat out wrong without the fix):
- [ ] Update `docs/whatsnew.rst`
- [ ] Update the released version on the site to this version, adjusting the tag and note that this is a critical bugfix which corrects the X.Y release
- [ ] Update the `omnia-md/conda-recipies` repo to point at the corrected version
