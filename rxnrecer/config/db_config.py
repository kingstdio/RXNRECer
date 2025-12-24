"""
Database-based configuration management for RXNRECer
Supports querying configuration items from database and associating with feature lists
"""
import os
import json
import logging
from typing import Dict, List, Optional, Any
from pathlib import Path
from rxnrecer.utils.db_utils import DatabaseConnection, get_db_connection_from_env

logger = logging.getLogger(__name__)


class DatabaseConfigManager:
    """Manage configuration items from database"""
    
    def __init__(self, db_connection: Optional[DatabaseConnection] = None, 
                 config_table: str = 'config_items',
                 feature_table: str = 'features',
                 feature_config_table: str = 'feature_configs'):
        """
        Initialize database configuration manager
        
        Args:
            db_connection: DatabaseConnection instance. If None, will try to create from env vars
            config_table: Name of the configuration items table
            feature_table: Name of the features table
            feature_config_table: Name of the feature-config association table
        """
        self.db = db_connection or get_db_connection_from_env()
        self.config_table = config_table
        self.feature_table = feature_table
        self.feature_config_table = feature_config_table
        
        # Cache for configuration items
        self._config_cache: Dict[str, Any] = {}
        self._feature_cache: Dict[str, List[str]] = {}
    
    def ensure_tables(self):
        """Ensure required tables exist in the database"""
        if self.db.db_type == 'sqlite':
            self._create_sqlite_tables()
        elif self.db.db_type == 'mysql':
            self._create_mysql_tables()
        elif self.db.db_type in ('postgresql', 'postgres'):
            self._create_postgresql_tables()
    
    def _create_sqlite_tables(self):
        """Create tables for SQLite"""
        with self.db.get_connection() as conn:
            cursor = conn.cursor()
            
            # Config items table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.config_table} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    config_key TEXT UNIQUE NOT NULL,
                    config_value TEXT,
                    config_type TEXT DEFAULT 'string',
                    description TEXT,
                    category TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Features table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.feature_table} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    feature_name TEXT UNIQUE NOT NULL,
                    feature_description TEXT,
                    enabled BOOLEAN DEFAULT 1,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Feature-Config association table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.feature_config_table} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    feature_id INTEGER NOT NULL,
                    config_key TEXT NOT NULL,
                    FOREIGN KEY (feature_id) REFERENCES {self.feature_table}(id) ON DELETE CASCADE,
                    UNIQUE(feature_id, config_key)
                )
            """)
            
            conn.commit()
            cursor.close()
    
    def _create_mysql_tables(self):
        """Create tables for MySQL"""
        with self.db.get_connection() as conn:
            cursor = conn.cursor()
            
            # Config items table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.config_table} (
                    id INT AUTO_INCREMENT PRIMARY KEY,
                    config_key VARCHAR(255) UNIQUE NOT NULL,
                    config_value TEXT,
                    config_type VARCHAR(50) DEFAULT 'string',
                    description TEXT,
                    category VARCHAR(100),
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
                    INDEX idx_key (config_key),
                    INDEX idx_category (category)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
            """)
            
            # Features table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.feature_table} (
                    id INT AUTO_INCREMENT PRIMARY KEY,
                    feature_name VARCHAR(255) UNIQUE NOT NULL,
                    feature_description TEXT,
                    enabled BOOLEAN DEFAULT TRUE,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    INDEX idx_enabled (enabled)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
            """)
            
            # Feature-Config association table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.feature_config_table} (
                    id INT AUTO_INCREMENT PRIMARY KEY,
                    feature_id INT NOT NULL,
                    config_key VARCHAR(255) NOT NULL,
                    FOREIGN KEY (feature_id) REFERENCES {self.feature_table}(id) ON DELETE CASCADE,
                    UNIQUE KEY unique_feature_config (feature_id, config_key),
                    INDEX idx_config_key (config_key)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4
            """)
            
            conn.commit()
            cursor.close()
    
    def _create_postgresql_tables(self):
        """Create tables for PostgreSQL"""
        with self.db.get_connection() as conn:
            cursor = conn.cursor()
            
            # Config items table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.config_table} (
                    id SERIAL PRIMARY KEY,
                    config_key VARCHAR(255) UNIQUE NOT NULL,
                    config_value TEXT,
                    config_type VARCHAR(50) DEFAULT 'string',
                    description TEXT,
                    category VARCHAR(100),
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_key ON {self.config_table}(config_key)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_category ON {self.config_table}(category)")
            
            # Features table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.feature_table} (
                    id SERIAL PRIMARY KEY,
                    feature_name VARCHAR(255) UNIQUE NOT NULL,
                    feature_description TEXT,
                    enabled BOOLEAN DEFAULT TRUE,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_enabled ON {self.feature_table}(enabled)")
            
            # Feature-Config association table
            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.feature_config_table} (
                    id SERIAL PRIMARY KEY,
                    feature_id INTEGER NOT NULL,
                    config_key VARCHAR(255) NOT NULL,
                    FOREIGN KEY (feature_id) REFERENCES {self.feature_table}(id) ON DELETE CASCADE,
                    UNIQUE(feature_id, config_key)
                )
            """)
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_config_key ON {self.feature_config_table}(config_key)")
            
            conn.commit()
            cursor.close()
    
    def get_config(self, config_key: str, default: Any = None) -> Any:
        """
        Get a configuration value by key
        
        Args:
            config_key: Configuration key
            default: Default value if key not found
            
        Returns:
            Configuration value (converted to appropriate type)
        """
        # Check cache first
        if config_key in self._config_cache:
            return self._config_cache[config_key]
        
        query = f"""
            SELECT config_value, config_type 
            FROM {self.config_table} 
            WHERE config_key = ?
        """
        
        if self.db.db_type != 'sqlite':
            query = query.replace('?', '%s')
        
        results = self.db.execute_query(query, (config_key,))
        
        if not results:
            return default
        
        value = results[0]['config_value']
        config_type = results[0]['config_type']
        
        # Convert value to appropriate type
        if value is None:
            return default
        
        if config_type == 'int':
            value = int(value)
        elif config_type == 'float':
            value = float(value)
        elif config_type == 'bool':
            value = value.lower() in ('true', '1', 'yes', 'on')
        elif config_type == 'json':
            value = json.loads(value)
        elif config_type == 'list':
            value = value.split(',') if value else []
        
        # Cache the value
        self._config_cache[config_key] = value
        return value
    
    def set_config(self, config_key: str, config_value: Any, 
                   config_type: str = 'string', description: str = None, 
                   category: str = None):
        """
        Set a configuration value
        
        Args:
            config_key: Configuration key
            config_value: Configuration value
            config_type: Type of value (string, int, float, bool, json, list)
            description: Description of the configuration
            category: Category of the configuration
        """
        # Convert value to string for storage
        if config_type == 'json':
            value_str = json.dumps(config_value)
        elif config_type == 'list':
            value_str = ','.join(str(v) for v in config_value)
        else:
            value_str = str(config_value)
        
        # Check if config exists
        check_query = f"SELECT id FROM {self.config_table} WHERE config_key = ?"
        if self.db.db_type != 'sqlite':
            check_query = check_query.replace('?', '%s')
        
        existing = self.db.execute_query(check_query, (config_key,))
        
        if existing:
            # Update existing
            update_query = f"""
                UPDATE {self.config_table} 
                SET config_value = ?, config_type = ?, description = ?, category = ?, updated_at = CURRENT_TIMESTAMP
                WHERE config_key = ?
            """
            if self.db.db_type != 'sqlite':
                update_query = update_query.replace('?', '%s')
            
            self.db.execute_update(update_query, (value_str, config_type, description, category, config_key))
        else:
            # Insert new
            insert_query = f"""
                INSERT INTO {self.config_table} (config_key, config_value, config_type, description, category)
                VALUES (?, ?, ?, ?, ?)
            """
            if self.db.db_type != 'sqlite':
                insert_query = insert_query.replace('?', '%s')
            
            self.db.execute_update(insert_query, (config_key, value_str, config_type, description, category))
        
        # Update cache
        self._config_cache[config_key] = config_value
    
    def get_all_configs(self, category: Optional[str] = None) -> Dict[str, Any]:
        """
        Get all configuration items
        
        Args:
            category: Optional category filter
            
        Returns:
            Dictionary of config_key -> config_value
        """
        if category:
            query = f"SELECT config_key, config_value, config_type FROM {self.config_table} WHERE category = ?"
            if self.db.db_type != 'sqlite':
                query = query.replace('?', '%s')
            results = self.db.execute_query(query, (category,))
        else:
            query = f"SELECT config_key, config_value, config_type FROM {self.config_table}"
            results = self.db.execute_query(query)
        
        configs = {}
        for row in results:
            key = row['config_key']
            value = row['config_value']
            config_type = row['config_type']
            
            # Convert value
            if config_type == 'int':
                value = int(value)
            elif config_type == 'float':
                value = float(value)
            elif config_type == 'bool':
                value = value.lower() in ('true', '1', 'yes', 'on')
            elif config_type == 'json':
                value = json.loads(value)
            elif config_type == 'list':
                value = value.split(',') if value else []
            
            configs[key] = value
            self._config_cache[key] = value
        
        return configs
    
    def add_feature(self, feature_name: str, feature_description: str = None, enabled: bool = True) -> int:
        """
        Add a feature to the feature list
        
        Args:
            feature_name: Name of the feature
            feature_description: Description of the feature
            enabled: Whether the feature is enabled
            
        Returns:
            Feature ID
        """
        # Check if feature exists
        check_query = f"SELECT id FROM {self.feature_table} WHERE feature_name = ?"
        if self.db.db_type != 'sqlite':
            check_query = check_query.replace('?', '%s')
        
        existing = self.db.execute_query(check_query, (feature_name,))
        
        if existing:
            # Update existing
            update_query = f"""
                UPDATE {self.feature_table} 
                SET feature_description = ?, enabled = ?
                WHERE feature_name = ?
            """
            if self.db.db_type != 'sqlite':
                update_query = update_query.replace('?', '%s')
            
            self.db.execute_update(update_query, (feature_description, enabled, feature_name))
            return existing[0]['id']
        else:
            # Insert new
            insert_query = f"""
                INSERT INTO {self.feature_table} (feature_name, feature_description, enabled)
                VALUES (?, ?, ?)
            """
            if self.db.db_type != 'sqlite':
                insert_query = insert_query.replace('?', '%s')
            
            self.db.execute_update(insert_query, (feature_name, feature_description, enabled))
            
            # Get the inserted ID
            result = self.db.execute_query(check_query, (feature_name,))
            return result[0]['id']
    
    def associate_config_with_feature(self, feature_name: str, config_key: str):
        """
        Associate a configuration item with a feature
        
        Args:
            feature_name: Name of the feature
            config_key: Configuration key to associate
        """
        # Get feature ID
        feature_query = f"SELECT id FROM {self.feature_table} WHERE feature_name = ?"
        if self.db.db_type != 'sqlite':
            feature_query = feature_query.replace('?', '%s')
        
        features = self.db.execute_query(feature_query, (feature_name,))
        if not features:
            raise ValueError(f"Feature '{feature_name}' not found")
        
        feature_id = features[0]['id']
        
        # Check if association exists
        check_query = f"SELECT id FROM {self.feature_config_table} WHERE feature_id = ? AND config_key = ?"
        if self.db.db_type != 'sqlite':
            check_query = check_query.replace('?', '%s')
        
        existing = self.db.execute_query(check_query, (feature_id, config_key))
        
        if not existing:
            # Create association
            insert_query = f"""
                INSERT INTO {self.feature_config_table} (feature_id, config_key)
                VALUES (?, ?)
            """
            if self.db.db_type != 'sqlite':
                insert_query = insert_query.replace('?', '%s')
            
            self.db.execute_update(insert_query, (feature_id, config_key))
    
    def get_feature_configs(self, feature_name: str) -> Dict[str, Any]:
        """
        Get all configuration items associated with a feature
        
        Args:
            feature_name: Name of the feature
            
        Returns:
            Dictionary of config_key -> config_value
        """
        query = f"""
            SELECT c.config_key, c.config_value, c.config_type
            FROM {self.config_table} c
            INNER JOIN {self.feature_config_table} fc ON c.config_key = fc.config_key
            INNER JOIN {self.feature_table} f ON fc.feature_id = f.id
            WHERE f.feature_name = ?
        """
        if self.db.db_type != 'sqlite':
            query = query.replace('?', '%s')
        
        results = self.db.execute_query(query, (feature_name,))
        
        configs = {}
        for row in results:
            key = row['config_key']
            value = row['config_value']
            config_type = row['config_type']
            
            # Convert value
            if config_type == 'int':
                value = int(value)
            elif config_type == 'float':
                value = float(value)
            elif config_type == 'bool':
                value = value.lower() in ('true', '1', 'yes', 'on')
            elif config_type == 'json':
                value = json.loads(value)
            elif config_type == 'list':
                value = value.split(',') if value else []
            
            configs[key] = value
        
        return configs
    
    def get_all_features(self, enabled_only: bool = False) -> List[Dict[str, Any]]:
        """
        Get all features
        
        Args:
            enabled_only: Only return enabled features
            
        Returns:
            List of feature dictionaries
        """
        if enabled_only:
            query = f"SELECT * FROM {self.feature_table} WHERE enabled = ? ORDER BY feature_name"
            if self.db.db_type != 'sqlite':
                query = query.replace('?', '%s')
            results = self.db.execute_query(query, (True,))
        else:
            query = f"SELECT * FROM {self.feature_table} ORDER BY feature_name"
            results = self.db.execute_query(query)
        
        return results
    
    def clear_cache(self):
        """Clear configuration cache"""
        self._config_cache.clear()
        self._feature_cache.clear()

