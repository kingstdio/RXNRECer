# RXNRECer Deployment Guide

## 🚀 Setting Release Branch as Default

This guide explains how to set the `release` branch as the default branch on GitHub and complete the release process.

### 📋 Prerequisites

- ✅ GitHub repository: https://github.com/kingstdio/RXNRECer
- ✅ Release branch created and pushed
- ✅ Release tag v1.0.0 created
- ✅ All necessary files committed and pushed

### 🔧 Step 1: Set Release Branch as Default

1. **Go to GitHub Repository Settings**
   - Navigate to: https://github.com/kingstdio/RXNRECer
   - Click on "Settings" tab

2. **Change Default Branch**
   - In the left sidebar, click "Branches"
   - Under "Default branch", click the dropdown
   - Select `release` from the list
   - Click "Update" to confirm
   - Click "I understand, update the default branch"

### 🔄 Step 2: Update Branch Protection (Optional but Recommended)

1. **Go to Branch Protection Rules**
   - In Settings → Branches, click "Add rule"
   - Branch name pattern: `release`

2. **Configure Protection Rules**
   - ✅ Require a pull request before merging
   - ✅ Require status checks to pass before merging
   - ✅ Require branches to be up to date before merging
   - ✅ Include administrators
   - ✅ Restrict pushes that create files
   - ✅ Restrict deletions

### 🏷️ Step 3: Create GitHub Release

1. **Go to Releases**
   - Click on "Releases" in the right sidebar
   - Click "Create a new release"

2. **Configure Release**
   - **Tag version**: `v1.0.0`
   - **Target**: `release` branch
   - **Release title**: `RXNRECer v1.0.0 - Production Release`
   - **Description**: Copy content from `RELEASE_NOTES.md`

3. **Upload Assets**
   - Upload the built package files (if available)
   - Or let GitHub Actions handle this automatically

4. **Publish Release**
   - Click "Publish release"

### 🔄 Step 4: Update Development Workflow

1. **Create Feature Branches from Release**
   ```bash
   git checkout release
   git pull origin release
   git checkout -b feature/new-feature
   ```

2. **Submit Pull Requests to Release Branch**
   - All new features should target the `release` branch
   - Use the `main` branch for experimental features only

### 📊 Step 5: Monitor and Maintain

1. **Check GitHub Actions**
   - Verify that CI/CD pipeline is working
   - Monitor test results and builds

2. **Update Documentation**
   - Keep README.md up to date
   - Update RELEASE_NOTES.md for new versions

3. **Version Management**
   - Create new tags for each release: `v1.1.0`, `v1.2.0`, etc.
   - Update setup.py version numbers accordingly

### 🚨 Important Notes

- **Never force push to the `release` branch** - it's now the default
- **Always create feature branches** from `release` for new development
- **Use pull requests** to merge changes into `release`
- **Test thoroughly** before merging to `release`

### 🔗 Useful Links

- **Repository**: https://github.com/kingstdio/RXNRECer
- **Issues**: https://github.com/kingstdio/RXNRECer/issues
- **Actions**: https://github.com/kingstdio/RXNRECer/actions
- **Releases**: https://github.com/kingstdio/RXNRECer/releases

### 📞 Support

If you encounter any issues during deployment:

1. Check the GitHub Actions logs
2. Review the error messages
3. Contact: zhenkun.shi@tib.cas.cn
4. Create an issue on GitHub

---

**Congratulations! 🎉** You've successfully set up RXNRECer with a production-ready release branch as the default.
