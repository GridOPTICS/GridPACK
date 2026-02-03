# CI/CD Pipeline Documentation

## Overview

GridPACK uses GitHub Actions to build multi-architecture Docker images that are automatically published to Docker Hub at `pnnl/gridpack`.

## Pipeline Architecture

The pipeline uses native ARM64 and AMD64 runners for optimal build performance:

```
prepare-metadata → build-native (amd64 + arm64) → merge-manifest → test-native (amd64 + arm64)
                        ↓ parallel ↓                       ↓              ↓ parallel ↓
```

### Jobs

1. **prepare-metadata**: Generates Docker tags and labels
2. **build-native**: Builds images natively on AMD64 and ARM64 runners in parallel
3. **merge-manifest**: Stitches architecture-specific images into a multi-arch manifest
4. **test-native**: Runs ctest smoke tests on both architectures in parallel

### Performance

- **Native ARM64 builds**: 3-10x faster than QEMU emulation
- **Parallel execution**: Both architectures build simultaneously
- **Total build time**: ~20-35 minutes (vs 45-90 minutes with QEMU)

## Triggers

### Automatic: Release Published

When a GitHub release is created, the pipeline automatically builds and publishes:

```bash
gh release create v3.3.0 --title "Release 3.3.0" --notes "Release notes"
```

Generated tags:
- `pnnl/gridpack:3.3.0`
- `pnnl/gridpack:3.3`
- `pnnl/gridpack:latest`

### Manual: Workflow Dispatch

Trigger builds manually from the GitHub Actions UI or CLI:

**Via GitHub UI:**
1. Go to: Actions → Docker Multi-Architecture Build and Publish
2. Click "Run workflow"
3. Configure options:
   - **tag**: Custom Docker tag (default: `dev`)
   - **skip_tests**: Skip ctest smoke tests (default: `false`)
   - **dockerfile**: Dockerfile path (default: `./Dockerfile`)
4. Click "Run workflow"

**Via GitHub CLI:**

```bash
# Basic run with defaults
gh workflow run docker-build.yml

# With custom tag
gh workflow run docker-build.yml -f tag=test-build

# Skip tests for faster iteration
gh workflow run docker-build.yml -f skip_tests=true

# Use custom Dockerfile
gh workflow run docker-build.yml -f dockerfile=./Dockerfile.custom
```

## Docker Image

### Repository
- **Location**: `pnnl/gridpack`
- **Registry**: Docker Hub

### Supported Architectures
- `linux/amd64` (Intel/AMD)
- `linux/arm64` (ARM64/Apple Silicon/AWS Graviton)

### Usage

See [Docker Usage Guide](./DOCKER.md) for detailed examples.

## Testing

The pipeline runs ctest smoke tests on both architectures after building. Key points:

- **Non-blocking**: Test failures do NOT prevent image publication (`continue-on-error: true`)
- **Expected pass rate**: ~94%
- **Can be skipped**: Use `skip_tests=true` for faster iteration

## Secrets Required

The pipeline requires two GitHub repository secrets:

| Secret | Description |
|--------|-------------|
| `DOCKERHUB_USERNAME` | Docker Hub username for authentication |
| `DOCKERHUB_TOKEN` | Docker Hub access token (not password!) |

### Setting Up Secrets

1. **Generate Docker Hub Access Token:**
   - Go to: https://hub.docker.com → Account Settings → Security
   - Create token: Name: `gridpack-ci`, Permissions: Read, Write, Delete
   - Copy the token (shown only once)

2. **Add to GitHub:**
   - Go to: Repository Settings → Secrets and variables → Actions
   - Add both secrets with the values above

## Monitoring

### View Workflow Runs

```bash
# List recent runs
gh run list --workflow=docker-build.yml --limit 5

# Watch a specific run
gh run watch <run-id>

# View detailed logs
gh run view <run-id> --log
```

### Inspect Published Image

```bash
# View manifest with both architectures
docker buildx imagetools inspect pnnl/gridpack:latest

# Expected output shows:
# - linux/amd64 manifest
# - linux/arm64 manifest
```

## Caching

The pipeline uses GitHub Actions cache to speed up builds:

- Separate cache per architecture (`build-amd64`, `build-arm64`)
- Automatically managed by GitHub
- Improves rebuild times significantly

## Troubleshooting

### Build Fails on One Architecture

The pipeline uses `fail-fast: false`, so if one architecture fails, the other continues. Check job logs for architecture-specific issues.

### Manifest Merge Fails

This usually indicates both architecture builds failed. Check that:
- Docker Hub credentials are correct
- Both build jobs completed successfully
- Digest artifacts were uploaded

### Authentication Issues

Verify secrets are configured correctly:
```bash
gh secret list
```

Should show:
- `DOCKERHUB_USERNAME`
- `DOCKERHUB_TOKEN`

## Technical Details

### Push-by-Digest Workflow

The pipeline uses a digest-based approach:

1. Each architecture builds and pushes: `pnnl/gridpack@sha256:abc123...`
2. These are content-addressed (by digest), not tagged yet
3. Manifest merge creates tags that reference both digests
4. Docker automatically selects the right image based on client platform

### Why Native Runners?

- **Performance**: 3-10x faster ARM64 builds without QEMU
- **Reliability**: No emulation issues or crashes
- **Testing**: Tests run on actual target platforms
- **Cost**: Same GitHub Actions pricing as AMD64

## Workflow File

Location: `.github/workflows/docker-build.yml`

The workflow definition must be in the default branch to be accessible. Changes to the workflow require merging to the default branch.

## Future Enhancements

Potential improvements:
- Additional architectures (RISC-V when runners available)
- Multi-registry support (GitHub Container Registry)
- Build time metrics/reporting
- Automated security scanning per architecture
