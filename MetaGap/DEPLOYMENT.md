# MetaGap Production Deployment Guide

This guide provides step-by-step instructions for deploying the MetaGap Django application using Nginx and Gunicorn on a dnf-based Linux system (Fedora, RHEL, CentOS, AlmaLinux, Rocky Linux).

## Prerequisites

- Root or sudo access to your server
- Domain name pointing to your server (optional but recommended)
- MetaGap application code in `/workspace/MetaGap`

## Quick Deploy

Run the automated deployment script:

```bash
chmod +x /workspace/MetaGap/deploy.sh
sudo /workspace/MetaGap/deploy.sh
```

## Manual Step-by-Step Setup

### Step 1: Install Required Packages

```bash
sudo dnf update -y
sudo dnf install -y python3-pip python3-devel gcc nginx
pip3 install gunicorn
```

### Step 2: Create Application User

```bash
sudo useradd -r -s /bin/false -d /var/www/metagap metagap
```

### Step 3: Set Up Application Directory

```bash
sudo mkdir -p /var/www/metagap
sudo cp -r /workspace/MetaGap/* /var/www/metagap/
sudo chown -R metagap:www-data /var/www/metagap
sudo chmod -R 755 /var/www/metagap
```

### Step 4: Install Python Dependencies

```bash
cd /var/www/metagap
pip3 install -r requirements.txt
```

### Step 5: Collect Static Files

```bash
cd /var/www/metagap
python3 manage.py collectstatic --noinput
```

### Step 6: Configure Gunicorn Systemd Service

The `gunicorn.service` file has been configured with:
- **User**: `metagap`
- **Group**: `www-data`
- **Working Directory**: `/var/www/metagap`
- **Socket**: `/run/gunicorn.sock`
- **WSGI Application**: `MetaGap.wsgi:application`

Install the service:

```bash
sudo cp /workspace/MetaGap/gunicorn.service /etc/systemd/system/gunicorn.service
sudo systemctl daemon-reload
sudo systemctl enable gunicorn
sudo systemctl start gunicorn
```

### Step 7: Configure Nginx

The `nginx.conf` file is configured with:
- **Static files**: Served directly from `/var/www/metagap/staticfiles/`
- **Media files**: Served directly from `/var/www/metagap/media/`
- **Proxy**: Passes requests to Gunicorn via Unix socket

**Important**: Update the `server_name` directive with your actual domain!

```bash
# Edit the nginx.conf file before copying
sudo nano /workspace/MetaGap/nginx.conf
# Change: server_name metagap.example.com www.metagap.example.com;
# To:     server_name your-domain.com www.your-domain.com;

# Then install the configuration
sudo cp /workspace/MetaGap/nginx.conf /etc/nginx/conf.d/metagap.conf
sudo nginx -t
sudo systemctl enable nginx
sudo systemctl start nginx
```

### Step 8: Configure Firewall

```bash
sudo firewall-cmd --permanent --add-service=http
sudo firewall-cmd --permanent --add-service=https
sudo firewall-cmd --reload
```

### Step 9: Verify Deployment

```bash
# Check service status
sudo systemctl status gunicorn
sudo systemctl status nginx

# Check socket permissions
ls -la /run/gunicorn.sock

# Test the application
curl http://localhost
```

## Configuration Files

### gunicorn.service

Located at: `/workspace/MetaGap/gunicorn.service`

Key settings:
- Binds to Unix socket: `/run/gunicorn.sock`
- Runs 3 workers
- Auto-restarts on failure
- Uses system-wide gunicorn installation

### nginx.conf

Located at: `/workspace/MetaGap/nginx.conf`

Key settings:
- Listens on port 80
- Serves static files with caching headers
- Proxies dynamic requests to Gunicorn
- Includes proper proxy headers

## Post-Deployment Tasks

1. **Update Domain Name**: Edit `/etc/nginx/conf.d/metagap.conf` and change `server_name` to your actual domain

2. **Configure SSL** (Recommended):
   ```bash
   sudo dnf install -y certbot python3-certbot-nginx
   sudo certbot --nginx -d your-domain.com
   ```

3. **Secure Environment Variables**:
   - Set `DEBUG=0` in `.env`
   - Generate a secure `SECRET_KEY`
   - Update `ALLOWED_HOSTS` with your production domains

4. **Set Up Database Backups**: Configure regular backups for your SQLite database or migrate to PostgreSQL for production

5. **Monitor Logs**:
   ```bash
   # Gunicorn logs
   journalctl -u gunicorn -f
   
   # Nginx logs
   sudo tail -f /var/log/nginx/access.log
   sudo tail -f /var/log/nginx/error.log
   ```

## Troubleshooting

### Gunicorn not starting
```bash
sudo systemctl status gunicorn
sudo journalctl -u gunicorn -n 50
```

### Nginx configuration errors
```bash
sudo nginx -t
sudo tail -f /var/log/nginx/error.log
```

### Permission issues
```bash
sudo chown -R metagap:www-data /var/www/metagap
sudo chmod 660 /run/gunicorn.sock
sudo chown metagap:www-data /run/gunicorn.sock
```

### Application errors
```bash
cd /var/www/metagap
python3 manage.py check
python3 manage.py runserver 0.0.0.0:8000  # For testing only
```

## File Structure

```
/var/www/metagap/
├── MetaGap/              # Django project settings
│   ├── settings.py
│   ├── urls.py
│   └── wsgi.py
├── app/                  # Main Django application
├── staticfiles/          # Collected static files
├── media/                # User-uploaded media files
├── db.sqlite3            # SQLite database
├── manage.py             # Django management script
├── requirements.txt      # Python dependencies
└── .env                  # Environment variables
```

## Support

For issues specific to MetaGap functionality, refer to the project documentation.
For deployment issues, check systemd and Nginx logs.
