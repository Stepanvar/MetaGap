const THEME_STORAGE_KEY = 'metagap-theme';
const LIGHT_THEME = 'light';
const DARK_THEME = 'dark';

const getStoredTheme = () => localStorage.getItem(THEME_STORAGE_KEY);

const storeTheme = (theme) => {
    localStorage.setItem(THEME_STORAGE_KEY, theme);
};

const getPreferredTheme = () => {
    const storedTheme = getStoredTheme();
    if (storedTheme) {
        return storedTheme;
    }

    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches
        ? DARK_THEME
        : LIGHT_THEME;
};

const updateToggleAppearance = (theme) => {
    const toggleButton = document.getElementById('themeToggle');
    if (!toggleButton) {
        return;
    }

    const lightIcon = toggleButton.querySelector('[data-icon-light]');
    const darkIcon = toggleButton.querySelector('[data-icon-dark]');
    const label = toggleButton.querySelector('[data-theme-label]');

    if (lightIcon && darkIcon) {
        if (theme === DARK_THEME) {
            lightIcon.classList.add('d-none');
            darkIcon.classList.remove('d-none');
        } else {
            lightIcon.classList.remove('d-none');
            darkIcon.classList.add('d-none');
        }
    }

    if (label) {
        label.textContent = theme === DARK_THEME ? 'Dark' : 'Light';
    }
};

const applyTheme = (theme) => {
    document.documentElement.setAttribute('data-bs-theme', theme);
    updateToggleAppearance(theme);
};

const initThemeToggle = () => {
    const toggleButton = document.getElementById('themeToggle');
    if (!toggleButton) {
        return;
    }

    toggleButton.addEventListener('click', () => {
        const currentTheme = document.documentElement.getAttribute('data-bs-theme');
        const newTheme = currentTheme === DARK_THEME ? LIGHT_THEME : DARK_THEME;
        applyTheme(newTheme);
        storeTheme(newTheme);
    });
};

document.addEventListener('DOMContentLoaded', () => {
    const preferredTheme = getPreferredTheme();
    applyTheme(preferredTheme);
    initThemeToggle();

    if (window.matchMedia) {
        const mediaQuery = window.matchMedia('(prefers-color-scheme: dark)');
        mediaQuery.addEventListener('change', (event) => {
            const storedTheme = getStoredTheme();
            if (storedTheme) {
                return;
            }

            const newTheme = event.matches ? DARK_THEME : LIGHT_THEME;
            applyTheme(newTheme);
        });
    }
});
