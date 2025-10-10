# MetaGaP interface tokens & layout guide

This document captures the design primitives introduced in `styles.css`. Reference it when building new screens to stay aligned with the shared system.

## Color & sizing tokens

`styles.css` exposes the following CSS variables at `:root` for light/dark themes. Use them instead of hard-coded values:

| Token | Purpose |
| --- | --- |
| `--color-primary`, `--color-secondary` | Brand accent colors sourced from the Django site settings. |
| `--color-primary-contrast`, `--color-secondary-contrast` | Legible text color that sits on top of the brand accents. |
| `--header-height` | Fixed header height, used to offset body padding. |
| `--app-container-max-width` | Maximum width for page content (`app-container`). |
| `--app-page-gap`, `--app-stack-gap`, `--app-section-gap` | Consistent spacing tokens for page-level stacks, nested sections, and intra-section rhythm. |
| `--app-surface-*` | Background, border, radius, and shadow values for elevated cards and panels. |

Theme-specific tokens (e.g., `--app-body-bg`, `--app-navbar-bg`) are set on `[data-bs-theme='light']`/`[data-bs-theme='dark']`. Switching themes automatically updates the shell, surfaces, breadcrumbs, and tables. Avoid redefining these colors directly in templates.

## Structural primitives

* **`.app-container`** – Responsive wrapper used in the navbar, breadcrumbs, main content, and footer. Provides horizontal padding and caps content width using `--app-container-max-width`.
* **`.app-main`** – Vertical padding wrapper for the main element. Automatically offsets the fixed header.
* **`.app-content-stack`** – Flex column that spaces page header, surfaces, and sections using `--app-stack-gap`.
* **`.app-page-header`**, `.app-page-title`, `.app-page-subtitle` – Standardised page heading block. Prefer overriding the `page_header` block in templates instead of duplicating heading markup.
* **`.app-section-stack`** and **`.app-section`** – Stack and section wrappers for content inside a surface. Sections automatically add borders between each other for visual hierarchy.

## Component utilities

* **`.app-card`** – Elevated card treatment that aligns with the translucent surfaces in light and dark mode.
* **`.app-card-grid`** – CSS grid helper for auto-fitting cards across breakpoints.
* **`.app-icon-circle`** – Circular icon badge with brand-aware background modifiers (`--primary`, `--success`).
* **`.app-definition-grid`** – Grid-based definition list for key-value account metadata.
* **`.app-list-group`** – Refined list-group styling for action-heavy rows; pairs with responsive button group rules.
* **`.app-form-wrapper`** – Wrap Django form partials to remove default margins and align vertical spacing.
* **`.app-table-wrapper`** – Surface wrapper for tables to maintain consistent radius, border, and box-shadow.
* **`.app-filter-card`** – Variant of `.app-card` used in search filters.

## Layout guidelines

1. **Breadcrumbs** – Prefer extending the `breadcrumb_items` block to append crumbs. If you must override the block, reuse the `.app-breadcrumbs` structure so spacing and blur stay consistent.
2. **Surfaces** – Place the main page content inside the default `app-surface` unless a page explicitly needs full-bleed elements. Introduce additional surfaces using `.app-section` to keep spacing consistent.
3. **Forms** – Wrap `partials/form.html` with `.app-form-wrapper` and lean on Bootstrap grid utilities for column layout. Avoid adding manual `<hr>` separators; rely on section stacking instead.
4. **Tables** – Enclose data tables inside `.app-table-wrapper` to inherit light/dark backgrounds and rounded corners. Pair with `.table-responsive` for mobile behaviour.
5. **Navigation** – Keep new header actions inside `.app-navbar-actions`. The CSS already handles wrapping and toggler alignment on narrow viewports.

## Accessibility & contrast

* All primary/danger/secondary buttons derive their contrast from the brand tokens, ensuring sufficient contrast in both themes.
* Cards, breadcrumbs, and the navbar use theme-specific backgrounds with opacity rather than fixed colours to preserve readability over the body gradient.
* When introducing new badges or accent colours, sample them against both light and dark theme backgrounds and prefer using `color-mix` with the brand tokens.

## Template checklist

Before adding a new page:

1. Override `page_header` (and optionally `breadcrumb_items`) instead of recreating headings.
2. Wrap the main content in `.app-section-stack` and break out logical groups with `.app-section`.
3. Use `.app-card` for repeatable card items and `.app-card-grid` for responsive layouts.
4. Include `.app-table-wrapper` for table-heavy pages.
5. Reference these tokens to avoid drifting spacing or colours.

Following this guide keeps future views consistent with the MetaGaP design language.
