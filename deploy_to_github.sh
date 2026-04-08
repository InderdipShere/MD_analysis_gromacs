#!/bin/bash

# Deploy Qn Order Parameter Project to GitHub
# Usage: ./deploy_to_github.sh <github-url>

if [ -z "$1" ]; then
    echo "Usage: ./deploy_to_github.sh <github-url>"
    echo "Example: ./deploy_to_github.sh https://github.com/InderdipShere/MD_analysis_gromacs.git"
    exit 1
fi

GITHUB_URL=$1
PROJECT_DIR="/Users/is284326/exe"

cd "$PROJECT_DIR" || exit 1

echo "==================================================="
echo "Qn Steinhardt Order Parameter - GitHub Deployment"
echo "==================================================="
echo ""
echo "Repository URL: $GITHUB_URL"
echo "Local directory: $PROJECT_DIR"
echo ""

# Check if git is initialized
if [ ! -d ".git" ]; then
    echo "❌ Error: Git repository not initialized"
    echo "Please run: git init"
    exit 1
fi

# Check current branch
CURRENT_BRANCH=$(git branch --show-current)
echo "Current branch: $CURRENT_BRANCH"
echo ""

# Check repository status
echo "Checking repository status..."
if git diff-index --quiet HEAD --; then
    echo "✅ Working directory is clean"
else
    echo "⚠️  Warning: Uncommitted changes detected"
    git status
    echo ""
fi

echo "Repository info:"
echo "  Commits: $(git rev-list --count HEAD)"
echo "  Files tracked: $(git ls-files | wc -l)"
echo ""

# Add remote
echo "Setting up remote..."
if git remote get-url origin > /dev/null 2>&1; then
    echo "⚠️  Remote 'origin' already exists"
    echo "Current remote: $(git remote get-url origin)"
    read -p "Replace with new URL? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git remote remove origin
        git remote add origin "$GITHUB_URL"
        echo "✅ Remote updated"
    fi
else
    git remote add origin "$GITHUB_URL"
    echo "✅ Remote added"
fi

# Set main branch
echo "Setting main as default branch..."
git branch -M main
echo "✅ Main branch set"

# Push
echo ""
echo "Pushing to GitHub..."
git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "==================================================="
    echo "✅ DEPLOYMENT SUCCESSFUL!"
    echo "==================================================="
    echo ""
    echo "Your repository is now on GitHub:"
    echo "  URL: $GITHUB_URL"
    echo ""
    echo "Next steps:"
    echo "  1. Visit: $GITHUB_URL"
    echo "  2. Add collaborators if needed"
    echo "  3. Create releases/tags for versions"
    echo "  4. Enable GitHub Pages (optional)"
    echo "  5. Set up discussions (optional)"
    echo ""
else
    echo ""
    echo "❌ Push failed!"
    echo "Please check your GitHub credentials and URL"
    exit 1
fi
